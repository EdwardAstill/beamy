from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List, Set, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from beamy.beam1d.analysis import LoadedBeam

from .frame import Frame
from ..core.loads import (
    FrameLoadCase,
    PointForce, Moment, DistributedForce, LoadCase
)
from .member import Member
from .node import Node
from .solver import analyze_frame_geometry, solve_displacements
from .results import MemberResultsDirect


def _collect_fixed_dofs(frame: Frame, node_to_idx: Dict[str, int]) -> tuple[List[int], Set[str]]:
    """
    Collect fixed DOFs from node supports and member-end constraints.

    Member constraints are a 12-digit string:
      - first 6 -> constraints at start node DOFs
      - last 6  -> constraints at end node DOFs
    A '1' means the DOF is fixed (same meaning as Node.support digit).

    Returns:
      - fixed_dofs: sorted list of global DOF indices to constrain
      - constrained_nodes: set of node IDs that have at least one constrained DOF
    """
    fixed_dofs_set: Set[int] = set()
    constrained_nodes: Set[str] = set()

    for nid, node in frame.nodes.items():
        if node.support:
            for d in range(6):
                if node.support[d] == "1":
                    fixed_dofs_set.add(node_to_idx[nid] * 6 + d)
                    constrained_nodes.add(nid)

    for member in frame.members:
        if member.constraints:
            for d in range(6):
                if member.constraints[d] == "1":
                    fixed_dofs_set.add(node_to_idx[member.start_node_id] * 6 + d)
                    constrained_nodes.add(member.start_node_id)
                if member.constraints[6 + d] == "1":
                    fixed_dofs_set.add(node_to_idx[member.end_node_id] * 6 + d)
                    constrained_nodes.add(member.end_node_id)

    fixed_dofs = sorted(fixed_dofs_set)
    return fixed_dofs, constrained_nodes

def _auto_fix_rotations_for_truss_only_nodes(frame: Frame, node_to_idx: Dict[str, int], fixed_dofs: List[int]) -> List[int]:
    """
    Frame DOFs include rotations at each node. For truss/cable-only nodes, rotations are
    physically irrelevant and have no stiffness contribution, which can make K singular.

    To keep the linear system solvable, automatically fix RX/RY/RZ at nodes that are not
    connected to any beam element.
    """
    fixed_set = set(fixed_dofs)
    for nid, node in frame.nodes.items():
        has_beam = False
        for mid in node.connected_members:
            m = frame.get_member(mid)
            if m.element_type == "beam":
                has_beam = True
                break
        if not has_beam:
            base = node_to_idx[nid] * 6
            fixed_set.add(base + 3)
            fixed_set.add(base + 4)
            fixed_set.add(base + 5)
    return sorted(fixed_set)

@dataclass
class LoadedFrame:
    """A frame with loads applied, ready for analysis."""
    frame: Frame
    loads: FrameLoadCase
    
    _nodal_displacements: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _reactions: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(init=False, default_factory=dict)

    # Auto-splitting / bundling metadata
    original_frame: Frame = field(init=False)
    original_loads: FrameLoadCase = field(init=False)
    _member_bundle: Dict[str, List[str]] = field(init=False, default_factory=dict)  # parent -> ordered segment ids
    _member_segment_parent: Dict[str, str] = field(init=False, default_factory=dict)  # segment -> parent
    _member_segment_offset: Dict[str, float] = field(init=False, default_factory=dict)  # segment -> x0 along parent
    _member_nodes_along: Dict[str, List[Tuple[float, str]]] = field(init=False, default_factory=dict)  # parent -> [(x,node_id)]
    
    def __post_init__(self) -> None:
        # Preserve user-level model for member bundling and conversion.
        self.original_frame = self.frame
        self.original_loads = self.loads

        # Insert real nodes at member-point attachment points (loads/moments/supports) and
        # split members into solver segments.
        self.frame, self.loads = self._expand_model_for_attachment_points(self.original_frame, self.original_loads)

        self._validate_loads()
        self._adjust_full_length_loads()
        self._analyze_frame()

    @staticmethod
    def _merge_support(a: Optional[str], b: Optional[str]) -> Optional[str]:
        if not a:
            return b
        if not b:
            return a
        return "".join("1" if (ca == "1" or cb == "1") else "0" for ca, cb in zip(a, b))

    @staticmethod
    def _coord_key(p: np.ndarray, tol: float = 1e-9) -> tuple[float, float, float]:
        return (
            round(float(p[0]) / tol) * tol,
            round(float(p[1]) / tol) * tol,
            round(float(p[2]) / tol) * tol,
        )

    def _get_reference_member(self, member_id: str) -> Member:
        """Resolve a member for local-axis reference.

        After auto-splitting, user scripts may still reference original member IDs.
        Those are available on self.original_frame.
        """
        try:
            return self.frame.get_member(member_id)
        except KeyError:
            return self.original_frame.get_member(member_id)
    
    def _validate_loads(self) -> None:
        for nf in self.loads.nodal_forces:
            if nf.node_id not in self.frame.nodes:
                raise ValueError(f"Unknown node '{nf.node_id}'")
            if getattr(nf, "coords", "global") == "local":
                if not nf.reference_member_id:
                    raise ValueError(f"NodalForce at '{nf.node_id}' uses local coords but has no reference_member_id")
                try:
                    _ = self._get_reference_member(nf.reference_member_id)
                except KeyError:
                    raise ValueError(f"Unknown reference member '{nf.reference_member_id}' for nodal force at '{nf.node_id}'")
        for nm in self.loads.nodal_moments:
            if nm.node_id not in self.frame.nodes:
                raise ValueError(f"Unknown node '{nm.node_id}'")
            if getattr(nm, "coords", "global") == "local":
                if not nm.reference_member_id:
                    raise ValueError(f"NodalMoment at '{nm.node_id}' uses local coords but has no reference_member_id")
                try:
                    _ = self._get_reference_member(nm.reference_member_id)
                except KeyError:
                    raise ValueError(f"Unknown reference member '{nm.reference_member_id}' for nodal moment at '{nm.node_id}'")
        mids = [m.id for m in self.frame.members]
        for mpf in self.loads.member_point_forces:
            if mpf.member_id not in mids:
                raise ValueError(f"Unknown member '{mpf.member_id}'")
        for mpm in getattr(self.loads, "member_point_moments", []):
            if mpm.member_id not in mids:
                raise ValueError(f"Unknown member '{mpm.member_id}'")
        for mdf in self.loads.member_distributed_forces:
            if mdf.member_id not in mids:
                raise ValueError(f"Unknown member '{mdf.member_id}'")

    def _expand_model_for_attachment_points(self, frame: Frame, loads: FrameLoadCase) -> tuple[Frame, FrameLoadCase]:
        """Insert real nodes at member-point attachment locations and split members.

        - Creates nodes at: point forces, point moments, point supports, distributed-load endpoints.
        - Splits each original member into segment members for the solver.
        - Applies member-wide supports to every node on that member after splitting.
        - Rewrites loads: point forces/moments become nodal loads; distributed loads are split per segment.

        Returns:
            (expanded_frame, expanded_loads)
        """

        # Collect split positions per original member
        split_pos: Dict[str, Set[float]] = {}

        def add_pos(mid: str, x: float) -> None:
            split_pos.setdefault(mid, set()).add(float(x))

        for m in frame.members:
            add_pos(m.id, 0.0)
            add_pos(m.id, m.length)

        for mpf in loads.member_point_forces:
            m = frame.get_member(mpf.member_id)
            x = mpf.position * m.length if mpf.position_type == "relative" else mpf.position
            add_pos(m.id, x)

        for mpm in getattr(loads, "member_point_moments", []):
            m = frame.get_member(mpm.member_id)
            x = mpm.position * m.length if mpm.position_type == "relative" else mpm.position
            add_pos(m.id, x)

        for mps in getattr(loads, "member_point_supports", []):
            m = frame.get_member(mps.member_id)
            x = mps.position * m.length if mps.position_type == "relative" else mps.position
            add_pos(m.id, x)

        for mdf in loads.member_distributed_forces:
            m = frame.get_member(mdf.member_id)
            s = float(mdf.start_position)
            e = float(m.length if mdf.end_position == -1.0 else mdf.end_position)
            add_pos(m.id, s)
            add_pos(m.id, e)

        # Connectivity: split any member at other members' end nodes that lie on it.
        # This allows user scripts to define continuous members without manually splitting
        # them just to create connection nodes.
        tol = 1e-9
        node_positions = list(frame.node_positions.values())
        for m in frame.members:
            start = frame.get_node(m.start_node_id).position
            dir_vec = m.direction
            L = m.length
            for p in node_positions:
                v = p - start
                x = float(np.dot(v, dir_vec))
                if x < -1e-9 or x > L + 1e-9:
                    continue
                perp = v - x * dir_vec
                if float(np.linalg.norm(perp)) <= tol:
                    add_pos(m.id, min(max(x, 0.0), L))

        needs_split = any(len(v) > 2 for v in split_pos.values()) or bool(getattr(loads, "member_supports", [])) or bool(getattr(loads, "member_point_supports", []))
        if not needs_split:
            # Trivial mapping for conversion.
            self._member_bundle = {m.id: [m.id] for m in frame.members}
            self._member_segment_parent = {m.id: m.id for m in frame.members}
            self._member_segment_offset = {m.id: 0.0 for m in frame.members}
            self._member_nodes_along = {m.id: [(0.0, m.start_node_id), (m.length, m.end_node_id)] for m in frame.members}
            return frame, loads

        # Node pool: start with existing nodes (preserve IDs)
        nodes_by_id: Dict[str, Node] = {}
        coord_to_node_id: Dict[tuple[float, float, float], str] = {}
        for nid, n in frame.nodes.items():
            nodes_by_id[nid] = Node(id=nid, position=n.position.copy(), support=n.support)
            coord_to_node_id[self._coord_key(n.position)] = nid

        # Next node id index
        max_n = -1
        for nid in nodes_by_id.keys():
            if nid.startswith("N"):
                try:
                    max_n = max(max_n, int(nid[1:]))
                except ValueError:
                    pass
        next_node_idx = max_n + 1

        def get_or_create_node_id(pos: np.ndarray) -> str:
            nonlocal next_node_idx
            key = self._coord_key(pos)
            if key in coord_to_node_id:
                return coord_to_node_id[key]
            nid = f"N{next_node_idx}"
            next_node_idx += 1
            nodes_by_id[nid] = Node(id=nid, position=pos.copy(), support=None)
            coord_to_node_id[key] = nid
            return nid

        # Segment members + bundle metadata
        new_members: List[Member] = []
        used_member_ids: Set[str] = set(m.id for m in frame.members)
        self._member_bundle = {}
        self._member_segment_parent = {}
        self._member_segment_offset = {}
        self._member_nodes_along = {}

        for parent in frame.members:
            L = parent.length
            xs_raw = sorted(split_pos.get(parent.id, {0.0, L}))
            # Clamp and de-dup
            xs: List[float] = []
            for x in xs_raw:
                xc = min(max(float(x), 0.0), L)
                if not xs or abs(xc - xs[-1]) > 1e-9:
                    xs.append(xc)
            if abs(xs[-1] - L) > 1e-9:
                xs.append(L)

            start_node = frame.get_node(parent.start_node_id)
            dir_vec = parent.direction
            node_chain: List[Tuple[float, str]] = []
            for x in xs:
                pos = start_node.position + dir_vec * x
                node_id = get_or_create_node_id(pos)
                node_chain.append((x, node_id))
            self._member_nodes_along[parent.id] = node_chain

            seg_ids: List[str] = []
            for i in range(len(node_chain) - 1):
                x0, n0 = node_chain[i]
                x1, n1 = node_chain[i + 1]
                if abs(x1 - x0) < 1e-12:
                    continue

                seg_id_base = f"{parent.id}__{i}"
                seg_id = seg_id_base
                k = 1
                while seg_id in used_member_ids:
                    k += 1
                    seg_id = f"{seg_id_base}_{k}"
                used_member_ids.add(seg_id)

                releases = None
                if parent.releases:
                    s_rel = parent.releases[0:6] if i == 0 else "000000"
                    e_rel = parent.releases[6:12] if i == (len(node_chain) - 2) else "000000"
                    releases = s_rel + e_rel

                constraints = None
                if parent.constraints:
                    s_c = parent.constraints[0:6] if i == 0 else "000000"
                    e_c = parent.constraints[6:12] if i == (len(node_chain) - 2) else "000000"
                    constraints = s_c + e_c

                new_members.append(
                    Member(
                        id=seg_id,
                        start_node_id=n0,
                        end_node_id=n1,
                        section=parent.section,
                        material=parent.material,
                        orientation=parent.orientation.copy(),
                        element_type=parent.element_type,
                        releases=releases,
                        constraints=constraints,
                    )
                )
                seg_ids.append(seg_id)
                self._member_segment_parent[seg_id] = parent.id
                self._member_segment_offset[seg_id] = float(x0)

            if not seg_ids:
                # Degenerate: keep original
                seg_ids = [parent.id]
                self._member_segment_parent[parent.id] = parent.id
                self._member_segment_offset[parent.id] = 0.0

            self._member_bundle[parent.id] = seg_ids

        # Apply member-point supports
        for mps in getattr(loads, "member_point_supports", []):
            parent = frame.get_member(mps.member_id)
            x = mps.position * parent.length if mps.position_type == "relative" else mps.position
            chain = self._member_nodes_along[parent.id]
            nearest_x, nid = min(chain, key=lambda t: abs(t[0] - x))
            if abs(nearest_x - x) > 1e-6:
                raise ValueError(f"MemberPointSupport on {parent.id} at x={x} did not land on a split node")
            nodes_by_id[nid].support = self._merge_support(nodes_by_id[nid].support, mps.support)

        # Apply member-wide supports to all nodes on that member
        for ms in getattr(loads, "member_supports", []):
            if ms.member_id not in self._member_nodes_along:
                raise ValueError(f"Unknown member '{ms.member_id}' for member support")
            for _x, nid in self._member_nodes_along[ms.member_id]:
                nodes_by_id[nid].support = self._merge_support(nodes_by_id[nid].support, ms.support)

        expanded_frame = Frame.from_nodes_and_members(list(nodes_by_id.values()), new_members, mpcs=frame.mpcs)

        # Rewrite loads
        expanded_loads = FrameLoadCase(name=loads.name)
        expanded_loads.nodal_forces = list(loads.nodal_forces)
        expanded_loads.nodal_moments = list(loads.nodal_moments)

        # Point loads/moments become nodal loads at the inserted node
        for mpf in loads.member_point_forces:
            parent = frame.get_member(mpf.member_id)
            x = mpf.position * parent.length if mpf.position_type == "relative" else mpf.position
            chain = self._member_nodes_along[parent.id]
            nearest_x, nid = min(chain, key=lambda t: abs(t[0] - x))
            if abs(nearest_x - x) > 1e-6:
                raise ValueError(f"MemberPointForce on {parent.id} at x={x} did not land on a split node")
            f_g = parent.transformation_matrix.T @ mpf.force if mpf.coords == "local" else mpf.force
            expanded_loads.add_nodal_force(nid, f_g, coords="global")

        for mpm in getattr(loads, "member_point_moments", []):
            parent = frame.get_member(mpm.member_id)
            x = mpm.position * parent.length if mpm.position_type == "relative" else mpm.position
            chain = self._member_nodes_along[parent.id]
            nearest_x, nid = min(chain, key=lambda t: abs(t[0] - x))
            if abs(nearest_x - x) > 1e-6:
                raise ValueError(f"MemberPointMoment on {parent.id} at x={x} did not land on a split node")
            m_g = parent.transformation_matrix.T @ mpm.moment if mpm.coords == "local" else mpm.moment
            expanded_loads.add_nodal_moment(nid, m_g, coords="global")

        # Distributed loads: split across segments (positions become segment-local)
        for mdf in loads.member_distributed_forces:
            parent = frame.get_member(mdf.member_id)
            L = parent.length
            s = float(mdf.start_position)
            e = float(L if mdf.end_position == -1.0 else mdf.end_position)
            if e < s:
                continue

            def f_at(x: float) -> np.ndarray:
                if abs(e - s) < 1e-12:
                    return mdf.start_force
                t = (x - s) / (e - s)
                return (1 - t) * mdf.start_force + t * mdf.end_force

            for seg_id in self._member_bundle[parent.id]:
                seg = expanded_frame.get_member(seg_id)
                x0 = float(self._member_segment_offset[seg_id])
                x1 = x0 + seg.length
                a = max(s, x0)
                b = min(e, x1)
                if b <= a + 1e-12:
                    continue
                expanded_loads.add_member_distributed_force(
                    member_id=seg_id,
                    start_position=a - x0,
                    end_position=b - x0,
                    start_force=f_at(a),
                    end_force=f_at(b),
                    coords=mdf.coords,
                )

        return expanded_frame, expanded_loads
    
    def _adjust_full_length_loads(self) -> None:
        for mdf in self.loads.member_distributed_forces:
            if mdf.end_position == -1.0:
                mdf.end_position = self.frame.get_member(mdf.member_id).length
    
    def _analyze_frame(self) -> None:
        node_ids = sorted(self.frame.nodes.keys())
        node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
        F_global = self._build_load_vector(node_to_idx, 6 * len(node_ids))
        
        fixed_dofs, constrained_nodes = _collect_fixed_dofs(self.frame, node_to_idx)
        fixed_dofs = _auto_fix_rotations_for_truss_only_nodes(self.frame, node_to_idx, fixed_dofs)

        cable_member_ids = [m.id for m in self.frame.members if m.element_type == "cable"]
        member_axial_scales: Dict[str, float] = {}
        for mid in cable_member_ids:
            member_axial_scales[mid] = 1.0

        # Tension-only cable iteration:
        # - Solve
        # - If a cable is in compression, slacken it (reduce axial stiffness)
        # - Repeat until cable set stops changing or max iterations reached
        slack_scale = 1e-9
        tol_n = 1e-6
        max_iter = 25

        d_global = None
        K_global = None
        m_mats = None

        for _it in range(max_iter):
            if cable_member_ids:
                K_global, m_mats = analyze_frame_geometry(self.frame, node_to_idx, member_axial_scales=member_axial_scales)
            else:
                K_global, m_mats = analyze_frame_geometry(self.frame, node_to_idx)

            d_global = solve_displacements(K_global, F_global, fixed_dofs)

            # Build per-node displacements for force recovery
            nodal_displacements: Dict[str, np.ndarray] = {}
            for nid in node_ids:
                nodal_displacements[nid] = d_global[node_to_idx[nid]*6 : node_to_idx[nid]*6 + 6]

            member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
            for member in self.frame.members:
                d_m_g = np.concatenate([nodal_displacements[member.start_node_id], nodal_displacements[member.end_node_id]])
                k_l, T = m_mats[member.id]
                f_m_l = k_l @ (T @ d_m_g)
                member_end_forces[member.id] = (f_m_l[0:6], f_m_l[6:12])

            if not cable_member_ids:
                self._nodal_displacements = nodal_displacements
                self._member_end_forces = member_end_forces
                break

            # Update cable activity
            changed = False
            for mid in cable_member_ids:
                start_f, _end_f = member_end_forces[mid]
                # Local axial force convention: tension => -startFx > 0
                n_axial = -float(start_f[0])
                desired = 1.0 if n_axial >= -tol_n else slack_scale
                if float(member_axial_scales[mid]) != float(desired):
                    member_axial_scales[mid] = float(desired)
                    changed = True

            if not changed:
                self._nodal_displacements = nodal_displacements
                self._member_end_forces = member_end_forces
                break

        if d_global is None or K_global is None:
            raise RuntimeError("Frame analysis did not run")

        R_global = K_global @ d_global - F_global
        for nid in constrained_nodes:
            self._reactions[nid] = R_global[node_to_idx[nid]*6 : node_to_idx[nid]*6 + 6]

    def _build_load_vector(self, node_to_idx: Dict[str, int], n_dofs: int) -> np.ndarray:
        F = np.zeros(n_dofs)
        for nf in self.loads.nodal_forces:
            f_g = nf.force
            if getattr(nf, "coords", "global") == "local":
                m_ref = self._get_reference_member(nf.reference_member_id)
                f_g = m_ref.transformation_matrix.T @ nf.force
            F[node_to_idx[nf.node_id]*6 : node_to_idx[nf.node_id]*6+3] += f_g

        for nm in self.loads.nodal_moments:
            m_g = nm.moment
            if getattr(nm, "coords", "global") == "local":
                m_ref = self._get_reference_member(nm.reference_member_id)
                m_g = m_ref.transformation_matrix.T @ nm.moment
            F[node_to_idx[nm.node_id]*6+3 : node_to_idx[nm.node_id]*6+6] += m_g
        for mpf in self.loads.member_point_forces:
            self._add_member_point_force_to_vector(mpf, node_to_idx, F)
        for mpm in getattr(self.loads, "member_point_moments", []):
            self._add_member_point_moment_to_vector(mpm, node_to_idx, F)
        for mdf in self.loads.member_distributed_forces:
            self._add_member_distributed_force_to_vector(mdf, node_to_idx, F)
        return F

    def _add_member_point_force_to_vector(self, mpf, node_to_idx, F):
        m = self.frame.get_member(mpf.member_id)
        x = mpf.position * m.length if mpf.position_type == "relative" else mpf.position
        f_g = m.transformation_matrix.T @ mpf.force if mpf.coords == "local" else mpf.force
        N1, N2 = 1 - x/m.length, x/m.length
        F[node_to_idx[m.start_node_id]*6 : node_to_idx[m.start_node_id]*6+3] += N1 * f_g
        F[node_to_idx[m.end_node_id]*6 : node_to_idx[m.end_node_id]*6+3] += N2 * f_g

    def _add_member_point_moment_to_vector(self, mpm, node_to_idx, F):
        m = self.frame.get_member(mpm.member_id)
        x = mpm.position * m.length if mpm.position_type == "relative" else mpm.position
        m_g = m.transformation_matrix.T @ mpm.moment if mpm.coords == "local" else mpm.moment
        N1, N2 = 1 - x/m.length, x/m.length
        F[node_to_idx[m.start_node_id]*6+3 : node_to_idx[m.start_node_id]*6+6] += N1 * m_g
        F[node_to_idx[m.end_node_id]*6+3 : node_to_idx[m.end_node_id]*6+6] += N2 * m_g

    def _add_member_distributed_force_to_vector(self, mdf, node_to_idx, F):
        m = self.frame.get_member(mdf.member_id)
        n_pts = 11
        xs = np.linspace(mdf.start_position, mdf.end_position, n_pts)
        dx = (mdf.end_position - mdf.start_position) / (n_pts - 1) if n_pts > 1 else 0
        for i, x in enumerate(xs):
            t = (x - mdf.start_position) / (mdf.end_position - mdf.start_position) if mdf.end_position > mdf.start_position else 0
            f_l = (1 - t) * mdf.start_force + t * mdf.end_force
            w = 0.5 if (i == 0 or i == n_pts - 1) else 1.0
            f_g = m.transformation_matrix.T @ (f_l * dx * w) if mdf.coords == "local" else (f_l * dx * w)
            N1, N2 = 1 - x/m.length, x/m.length
            F[node_to_idx[m.start_node_id]*6 : node_to_idx[m.start_node_id]*6+3] += N1 * f_g
            F[node_to_idx[m.end_node_id]*6 : node_to_idx[m.end_node_id]*6+3] += N2 * f_g

    @property
    def nodal_displacements(self) -> Dict[str, np.ndarray]: return self._nodal_displacements
    @property
    def reactions(self) -> Dict[str, np.ndarray]: return self._reactions
    @property
    def member_end_forces(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]: return self._member_end_forces

    def to_loaded_beams(self) -> Dict[str, 'LoadedBeam']:
        from ..beam1d.beam import Beam1D
        from ..core.support import Support
        from ..beam1d.analysis import LoadedBeam

        def _as_global(vec: np.ndarray, coords: str, reference_member_id: Optional[str]) -> np.ndarray:
            if coords != "local":
                return vec
            m_ref = self._get_reference_member(reference_member_id)
            return m_ref.transformation_matrix.T @ vec

        # If we didn't split, keep the legacy behavior.
        has_bundles = bool(self._member_bundle) and any(len(v) > 1 for v in self._member_bundle.values())
        if not has_bundles:
            loaded_beams: Dict[str, LoadedBeam] = {}
            for member in self.frame.members:
                start_f, end_f = self._member_end_forces[member.id]
                beam = Beam1D(L=member.length, material=member.material, section=member.section, supports=[Support(x=0.0, type="111111")])
                lc = LoadCase(name=f"Member {member.id} loads")

                if member.element_type != "beam":
                    lc.add_point_force(PointForce(point=np.array([member.length, 0, 0]), force=np.array([-end_f[0], 0, 0])))
                else:
                    lc.add_point_force(PointForce(point=np.array([member.length, 0, 0]), force=-end_f[0:3]))
                    if any(abs(end_f[3:6]) > 1e-10):
                        lc.add_moment(Moment(x=member.length, moment=-end_f[3:6]))
                    if any(abs(start_f[3:6]) > 1e-10):
                        lc.add_moment(Moment(x=0.0, moment=start_f[3:6]))

                for mpf in self.loads.member_point_forces:
                    if mpf.member_id == member.id:
                        x = mpf.position * member.length if mpf.position_type == "relative" else mpf.position
                        f_l = member.transformation_matrix @ mpf.force if mpf.coords == "global" else mpf.force
                        lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=f_l))

                for mpm in getattr(self.loads, "member_point_moments", []):
                    if mpm.member_id == member.id:
                        x = mpm.position * member.length if mpm.position_type == "relative" else mpm.position
                        m_l = member.transformation_matrix @ mpm.moment if mpm.coords == "global" else mpm.moment
                        lc.add_moment(Moment(x=x, moment=m_l))

                for mdf in self.loads.member_distributed_forces:
                    if mdf.member_id == member.id:
                        s_f_l = member.transformation_matrix @ mdf.start_force if mdf.coords == "global" else mdf.start_force
                        e_f_l = member.transformation_matrix @ mdf.end_force if mdf.coords == "global" else mdf.end_force
                        lc.add_distributed_force(
                            DistributedForce(
                                np.array([mdf.start_position, 0, 0]),
                                np.array([mdf.end_position, 0, 0]),
                                s_f_l,
                                e_f_l,
                            )
                        )

                loaded_beams[member.id] = LoadedBeam(beam, lc)
            return loaded_beams

        # Bundled path: one LoadedBeam per original member.
        loaded_beams: Dict[str, LoadedBeam] = {}
        for parent in self.original_frame.members:
            if parent.id not in self._member_bundle:
                continue
            seg_ids = self._member_bundle[parent.id]
            if not seg_ids:
                continue

            first_seg_id = seg_ids[0]
            last_seg_id = seg_ids[-1]
            start_f, _ = self._member_end_forces[first_seg_id]
            _, end_f = self._member_end_forces[last_seg_id]

            beam = Beam1D(L=parent.length, material=parent.material, section=parent.section, supports=[Support(x=0.0, type="111111")])
            lc = LoadCase(name=f"Member {parent.id} loads")

            if parent.element_type != "beam":
                lc.add_point_force(PointForce(point=np.array([parent.length, 0, 0]), force=np.array([-end_f[0], 0, 0])))
            else:
                lc.add_point_force(PointForce(point=np.array([parent.length, 0, 0]), force=-end_f[0:3]))
                if any(abs(end_f[3:6]) > 1e-10):
                    lc.add_moment(Moment(x=parent.length, moment=-end_f[3:6]))
                if any(abs(start_f[3:6]) > 1e-10):
                    lc.add_moment(Moment(x=0.0, moment=start_f[3:6]))

            chain = self._member_nodes_along.get(parent.id, [])
            node_to_x = {nid: float(x) for x, nid in chain}

            # Nodal loads on those nodes become beam point loads (in parent local coords)
            for nf in self.loads.nodal_forces:
                if nf.node_id not in node_to_x:
                    continue
                x = node_to_x[nf.node_id]
                f_g = _as_global(nf.force, getattr(nf, "coords", "global"), getattr(nf, "reference_member_id", None))
                f_l = parent.transformation_matrix @ f_g
                lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=f_l))

            for nm in self.loads.nodal_moments:
                if nm.node_id not in node_to_x:
                    continue
                x = node_to_x[nm.node_id]
                m_g = _as_global(nm.moment, getattr(nm, "coords", "global"), getattr(nm, "reference_member_id", None))
                m_l = parent.transformation_matrix @ m_g
                lc.add_moment(Moment(x=x, moment=m_l))

            # Distributed loads on segments become beam distributed loads (offset to parent x)
            for mdf in self.loads.member_distributed_forces:
                seg_id = mdf.member_id
                if seg_id not in self._member_segment_parent:
                    continue
                if self._member_segment_parent[seg_id] != parent.id:
                    continue
                x0 = float(self._member_segment_offset[seg_id])
                s_x = x0 + float(mdf.start_position)
                e_x = x0 + float(mdf.end_position)
                s_f_l = parent.transformation_matrix @ mdf.start_force if mdf.coords == "global" else mdf.start_force
                e_f_l = parent.transformation_matrix @ mdf.end_force if mdf.coords == "global" else mdf.end_force
                lc.add_distributed_force(DistributedForce(np.array([s_x, 0, 0]), np.array([e_x, 0, 0]), s_f_l, e_f_l))

            loaded_beams[parent.id] = LoadedBeam(beam, lc)

        return loaded_beams

    def get_member_results(self, member_id: str) -> MemberResultsDirect:
        """
        Get internal force results for a member using direct equilibrium computation.
        
        This correctly handles continuous members split at nodes by using the
        actual end forces from the global frame analysis.
        """
        member = self.frame.get_member(member_id)
        start_f, end_f = self._member_end_forces[member_id]
        
        # Collect distributed member loads in local coordinates
        member_loads = []
        for mdf in self.loads.member_distributed_forces:
            if mdf.member_id == member_id:
                # Convert to local coords if needed
                if mdf.coords == "global":
                    w_local = member.transformation_matrix @ mdf.start_force
                else:
                    w_local = mdf.start_force
                # For uniform loads, start_force == end_force
                member_loads.append((mdf.start_position, mdf.end_position, w_local))
        
        return MemberResultsDirect(
            member=member,
            start_f=start_f,
            end_f=end_f,
            member_loads=member_loads,
        )

    def to_loaded_beam(self, member_id: str) -> 'LoadedBeam':
        """
        Extract a single member (combining segments if split) as a LoadedBeam for detailed analysis.
        
        If the member was split during analysis, all segments are recombined into a single
        beam with their loads adjusted back to the original member's coordinate system.
        
        For nodes with fixed supports, those supports are preserved.
        For nodes without fixed supports, the reactions from the frame analysis
        are applied as boundary loads (opposite direction to balance the member).
        
        Args:
            member_id: The original (user-level) member ID to extract
            
        Returns:
            A LoadedBeam object with the member isolated and properly loaded
        """
        from ..beam1d.beam import Beam1D
        from ..core.support import Support
        from ..core.loads import LoadCase, PointForce, Moment, DistributedForce
        from ..beam1d.analysis import LoadedBeam
        
        # Resolve member: use original frame for user-level IDs
        parent = self.original_frame.get_member(member_id)
        
        # Get the segment IDs if this member was split
        seg_ids = self._member_bundle.get(member_id, [member_id])
        if not seg_ids:
            raise ValueError(f"No segments found for member {member_id}")
        
        first_seg_id = seg_ids[0]
        last_seg_id = seg_ids[-1]
        
        # Get end forces from the expanded frame analysis
        start_f, _ = self._member_end_forces[first_seg_id]
        _, end_f = self._member_end_forces[last_seg_id]
        
        # Build the isolated beam with start support
        beam = Beam1D(
            L=parent.length,
            material=parent.material,
            section=parent.section,
            supports=[Support(x=0.0, type="111111")]
        )
        lc = LoadCase(name=f"Member {member_id} extracted")
        
        # Get the node chain for this member to check for fixed supports
        chain = self._member_nodes_along.get(member_id, [])
        node_to_x = {nid: float(x) for x, nid in chain}
        
        # Find which nodes have fixed supports in the original frame
        fixed_node_ids = set()
        for nid, node in self.original_frame.nodes.items():
            if node.support and all(c == "1" for c in node.support):
                fixed_node_ids.add(nid)
        
        # Add loads from nodal forces/moments on the member's nodes
        for nf in self.loads.nodal_forces:
            if nf.node_id in node_to_x:
                x = node_to_x[nf.node_id]
                f_g = nf.force if getattr(nf, "coords", "global") == "global" else self._get_reference_member(getattr(nf, "reference_member_id", None)).transformation_matrix.T @ nf.force
                f_l = parent.transformation_matrix @ f_g
                lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=f_l))
        
        for nm in self.loads.nodal_moments:
            if nm.node_id in node_to_x:
                x = node_to_x[nm.node_id]
                m_g = nm.moment if getattr(nm, "coords", "global") == "global" else self._get_reference_member(getattr(nm, "reference_member_id", None)).transformation_matrix.T @ nm.moment
                m_l = parent.transformation_matrix @ m_g
                lc.add_moment(Moment(x=x, moment=m_l))
        
        # Add loads from member-level point forces/moments
        # Recombine loads across all segments back to original member coordinates
        for mpf in self.loads.member_point_forces:
            seg_id = mpf.member_id
            if seg_id not in self._member_segment_parent:
                continue
            if self._member_segment_parent[seg_id] != member_id:
                continue
            # Convert segment-local position to original member-local position
            x0 = float(self._member_segment_offset[seg_id])
            # mpf.position is relative to the segment; convert to original member coordinates
            x_in_segment = mpf.position * self.frame.get_member(seg_id).length if mpf.position_type == "relative" else mpf.position
            x = x0 + x_in_segment
            f_l = parent.transformation_matrix @ mpf.force if mpf.coords == "global" else mpf.force
            lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=f_l))
        
        for mpm in getattr(self.loads, "member_point_moments", []):
            seg_id = mpm.member_id
            if seg_id not in self._member_segment_parent:
                continue
            if self._member_segment_parent[seg_id] != member_id:
                continue
            # Convert segment-local position to original member-local position
            x0 = float(self._member_segment_offset[seg_id])
            x_in_segment = mpm.position * self.frame.get_member(seg_id).length if mpm.position_type == "relative" else mpm.position
            x = x0 + x_in_segment
            m_l = parent.transformation_matrix @ mpm.moment if mpm.coords == "global" else mpm.moment
            lc.add_moment(Moment(x=x, moment=m_l))
        
        # Add distributed member loads
        for mdf in self.loads.member_distributed_forces:
            seg_id = mdf.member_id
            if seg_id not in self._member_segment_parent:
                continue
            if self._member_segment_parent[seg_id] != member_id:
                continue
            # Convert segment-local positions to original member-local positions
            x0 = float(self._member_segment_offset[seg_id])
            s_in_segment = float(mdf.start_position)
            e_in_segment = float(mdf.end_position)
            s_x = x0 + s_in_segment
            e_x = x0 + e_in_segment
            s_f_l = parent.transformation_matrix @ mdf.start_force if mdf.coords == "global" else mdf.start_force
            e_f_l = parent.transformation_matrix @ mdf.end_force if mdf.coords == "global" else mdf.end_force
            lc.add_distributed_force(DistributedForce(np.array([s_x, 0, 0]), np.array([e_x, 0, 0]), s_f_l, e_f_l))
        
        # Apply reactions as boundary loads for non-fixed end nodes
        # Start node: already fixed (111111 support), don't add load
        # End node: apply negative end_f as load (reaction equilibrium)
        end_nid = parent.end_node_id
        
        if end_nid not in fixed_node_ids:
            # Apply end forces and moments as negative (equilibrium)
            lc.add_point_force(PointForce(point=np.array([parent.length, 0, 0]), force=-end_f[0:3]))
            if any(abs(end_f[3:6]) > 1e-10):
                lc.add_moment(Moment(x=parent.length, moment=-end_f[3:6]))
        
        # Check intermediate nodes for reactions and apply as loads
        for x, nid in chain:
            if x <= 1e-9 or x >= parent.length - 1e-9:
                continue  # Skip start/end nodes already handled
            if nid in self._reactions:
                rxn = self._reactions[nid]
                f_g = rxn[0:3]
                m_g = rxn[3:6]
                f_l = parent.transformation_matrix @ f_g
                m_l = parent.transformation_matrix @ m_g
                if any(abs(f_l) > 1e-10):
                    lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=-f_l))
                if any(abs(m_l) > 1e-10):
                    lc.add_moment(Moment(x=x, moment=-m_l))
        
        return LoadedBeam(beam, lc)

    def get_aisc_utilizations(self) -> Dict[str, float]:
        """
        Compute AISC Chapter F utilization ratios for all original members.

        Important:
            These utilizations are computed directly from the internal actions recovered from the
            *frame* analysis (end forces + member distributed loads). We do NOT re-solve each member
            as a separate 1D beam with artificial boundary conditions, because that can produce
            incorrect results for interior members (e.g. base cross members) and can break symmetry.

        Returns:
            Dictionary mapping original member ID to its governing utilization ratio (0-1 scale).
            Utilization = max(bending utilizations, shear utilization).
        """
        from ..beam1d.beam import Beam1D
        from ..core.support import Support
        from ..core.results import Result
        from ..checks import aisc_9

        utilizations: Dict[str, float] = {}

        class _AiscActionAdapter:
            """Minimal adapter to run AISC checks from precomputed internal actions."""

            def __init__(self, beam: Beam1D, my: Result, mz: Result, vz: Result):
                self.beam = beam
                self._my = my
                self._mz = mz
                self._vz = vz

            def bending(self, axis: str, points: int = 100):
                # aisc_9 expects:
                # - bending('y') -> strong-axis action paired with section.Sz
                # - bending('z') -> weak-axis action paired with section.Sy
                # In beamy frame results:
                # - bending_y is My
                # - bending_z is Mz
                if axis == "y":
                    return type("o", (), {"action": self._mz})()
                return type("o", (), {"action": self._my})()

            def shear(self, axis: str, points: int = 100):
                # aisc_9 currently uses shear('z') only.
                if axis == "z":
                    return type("o", (), {"action": self._vz})()
                return type("o", (), {"action": Result(self._vz._x, self._vz._values * 0.0)})()

        def _util_from_results(member: Member, my: Result, mz: Result, vz: Result) -> float:
            # Supports are only used by aisc_9 for unbraced-length segmentation.
            # Keep the existing conservative behavior: treat member as unbraced for full length.
            beam = Beam1D(
                L=member.length,
                material=member.material,
                section=member.section,
                supports=[Support(x=0.0, type="111111")],
            )
            adapter = _AiscActionAdapter(beam, my=my, mz=mz, vz=vz)
            return float(aisc_9.run(adapter, length_unit="m", force_unit="N").utilisation)

        def _stitch_parent_results(parent_id: str) -> tuple[Result, Result, Result]:
            parent = self.original_frame.get_member(parent_id)
            seg_ids = self._member_bundle.get(parent_id, [parent_id])

            # Sample the parent at a fixed grid so we can combine segments robustly.
            points = 801
            xs = np.linspace(0.0, parent.length, points)
            my = np.zeros_like(xs)
            mz = np.zeros_like(xs)
            vz = np.zeros_like(xs)

            offsets: List[Tuple[float, str, float]] = []
            for seg_id in seg_ids:
                seg = self.frame.get_member(seg_id)
                offsets.append((float(self._member_segment_offset.get(seg_id, 0.0)), seg_id, float(seg.length)))

            for i, x in enumerate(xs):
                for x0, seg_id, seg_L in offsets:
                    if x >= x0 - 1e-9 and x <= x0 + seg_L + 1e-9:
                        res = self.get_member_results(seg_id)
                        xl = float(x - x0)
                        my[i] = res.bending_y.action.at(xl)
                        mz[i] = res.bending_z.action.at(xl)
                        vz[i] = res.shear_z.action.at(xl)
                        break

            return Result(xs, my), Result(xs, mz), Result(xs, vz)

        # If we didn't split, compute per-solver-member (frame members are already original).
        has_bundles = bool(self._member_bundle) and any(len(v) > 1 for v in self._member_bundle.values())
        if not has_bundles:
            for member in self.frame.members:
                if member.element_type != "beam":
                    utilizations[member.id] = 0.0
                    continue
                try:
                    res = self.get_member_results(member.id)
                    util = _util_from_results(
                        member,
                        my=res.bending_y.action,
                        mz=res.bending_z.action,
                        vz=res.shear_z.action,
                    )
                    utilizations[member.id] = util
                except Exception:
                    utilizations[member.id] = 0.0
            return utilizations

        # Bundled path: one utilization per original member.
        for parent in self.original_frame.members:
            if parent.id not in self._member_bundle:
                continue
            if parent.element_type != "beam":
                utilizations[parent.id] = 0.0
                continue
            try:
                my, mz, vz = _stitch_parent_results(parent.id)
                utilizations[parent.id] = _util_from_results(parent, my=my, mz=mz, vz=vz)
            except Exception:
                utilizations[parent.id] = 0.0

        return utilizations

    def plot(self, **kwargs):
        from ..viz.frame_plots import plot_frame
        plot_frame(self, **kwargs)
    def plot_deflection(self, **kwargs):
        from ..viz.frame_plots import plot_deflection
        plot_deflection(self, **kwargs)
    def plot_von_mises(self, **kwargs):
        from ..viz.frame_plots import plot_von_mises
        plot_von_mises(self, **kwargs)
    def plot_results(self, **kwargs):
        from ..viz.frame_plots import plot_results
        plot_results(self, **kwargs)
    def plot_member_diagrams(self, member_id: str, **kwargs):
        from ..viz.frame_plots import plot_member_diagrams
        plot_member_diagrams(self, member_id, **kwargs)
    def plot_aisc_utilization(self, **kwargs):
        from ..viz.frame_plots import plot_aisc_utilization
        plot_aisc_utilization(self, **kwargs)
