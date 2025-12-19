from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List, Set, TYPE_CHECKING, Literal
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
from .solver import (
    assemble_global_stiffness,
    solve_displacements,
    recover_member_end_forces,
    ElementStiffnessScales,
    assemble_geometric_stiffness,
)
from ..core.math import build_transformation_matrix_12x12
from .results import MemberResultsDirect


AnalysisMethod = Literal[
    "GENERIC_LINEAR",
    "SECOND_ORDER_ELASTIC",
    "AISC360_DAM",
    "EC3_GLOBAL_IMPERFECTIONS",
    "AS4100_SECOND_ORDER",
]

ImperfectionModel = Literal["none", "notional_loads", "initial_sway"]
StiffnessRules = Literal["none", "aisc_dam", "ec3", "as4100"]
NotionalPSource = Literal["input_based", "reactions_based"]


@dataclass
class FrameAnalysisSettings:
    """Configuration for frame analysis.

    This is intentionally a first-class object so we can thread analysis choices
    (method, imperfections, stiffness rules, convergence controls) through the
    solver without retrofitting.
    """

    analysis_method: AnalysisMethod = "GENERIC_LINEAR"
    second_order: Optional[bool] = None

    imperfection_model: ImperfectionModel = "none"
    notional_factor: float = 0.002
    notional_p_source: NotionalPSource = "input_based"
    notional_axes: Tuple[Literal["x", "y"], ...] = ("x", "y")
    notional_signs: Tuple[float, float] = (1.0, 1.0)  # (x_sign, y_sign)

    stiffness_rules: StiffnessRules = "none"

    # Stiffness scaling (elastic modifiers).
    # If None, a default may be selected based on stiffness_rules.
    bending_stiffness_factor: Optional[float] = None  # scales Iy and Iz
    torsion_stiffness_factor: float = 1.0  # scales J
    axial_stiffness_factor: float = 1.0  # scales A

    # Iteration controls (used by nonlinear solvers; harmless for linear)
    max_iter: int = 25
    tol_u_rel: float = 1e-6
    tol_u_abs: float = 1e-12
    n_steps: int = 1
    relaxation_omega: float = 1.0

    # Existing nonlinear behavior
    cable_tension_only: bool = True

    def __post_init__(self) -> None:
        if self.second_order is None:
            self.second_order = self.analysis_method in (
                "SECOND_ORDER_ELASTIC",
                "AISC360_DAM",
                "EC3_GLOBAL_IMPERFECTIONS",
                "AS4100_SECOND_ORDER",
            )
        self.n_steps = max(int(self.n_steps), 1)
        self.max_iter = max(int(self.max_iter), 1)
        self.relaxation_omega = float(self.relaxation_omega)
        if not (0.0 < self.relaxation_omega <= 1.0):
            raise ValueError("relaxation_omega must be in (0, 1]")

        if self.notional_p_source not in ("input_based", "reactions_based"):
            raise ValueError("notional_p_source must be 'input_based' or 'reactions_based'")

        ax_ok = set(self.notional_axes).issubset({"x", "y"})
        if not ax_ok:
            raise ValueError("notional_axes must be a subset of ('x','y')")

        if len(self.notional_signs) != 2:
            raise ValueError("notional_signs must be a 2-tuple (x_sign, y_sign)")

        if self.bending_stiffness_factor is None:
            if self.stiffness_rules == "aisc_dam":
                self.bending_stiffness_factor = 0.8
            else:
                self.bending_stiffness_factor = 1.0

        self.bending_stiffness_factor = float(self.bending_stiffness_factor)
        self.torsion_stiffness_factor = float(self.torsion_stiffness_factor)
        self.axial_stiffness_factor = float(self.axial_stiffness_factor)

        if self.bending_stiffness_factor <= 0.0:
            raise ValueError("bending_stiffness_factor must be > 0")
        if self.torsion_stiffness_factor <= 0.0:
            raise ValueError("torsion_stiffness_factor must be > 0")
        if self.axial_stiffness_factor <= 0.0:
            raise ValueError("axial_stiffness_factor must be > 0")


@dataclass
class StabilizationReport:
    """Records numerical stabilization constraints added by the solver.

    These constraints are NOT physical supports and must not be interpreted as
    bracing/restraint for design checks.
    """

    kind: str
    added_dofs: List[int]
    node_ids: List[str]
    dof_names: List[str]

    @property
    def used(self) -> bool:
        return bool(self.added_dofs)


@dataclass
class FrameAnalysisResult:
    """Primary output container for frame analysis + checker-facing metadata."""

    nodal_displacements: Dict[str, np.ndarray]
    reactions: Dict[str, np.ndarray]
    member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]]

    # Metadata
    settings: FrameAnalysisSettings
    converged: bool
    iterations: int
    design_grade_ok: bool
    design_grade_notes: List[str]
    warnings: List[str]
    stabilization: StabilizationReport

    # Low-level info (useful for debugging)
    fixed_dofs_physical: List[int]
    fixed_dofs_total: List[int]


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


def _truss_only_rotation_stabilization(
    frame: Frame, node_to_idx: Dict[str, int], fixed_dofs: List[int]
) -> StabilizationReport:
    """Return stabilization constraints needed to avoid singular K.

    For truss/cable-only nodes, rotational DOFs can be unconstrained and have no stiffness.
    We fix RX/RY/RZ at those nodes as a numerical stabilizer.
    """
    fixed_set = set(fixed_dofs)
    added: List[int] = []
    nodes: List[str] = []

    for nid, node in frame.nodes.items():
        has_beam = False
        for mid in node.connected_members:
            m = frame.get_member(mid)
            if m.element_type == "beam":
                has_beam = True
                break
        if has_beam:
            continue

        base = node_to_idx[nid] * 6
        candidate = [base + 3, base + 4, base + 5]
        newly = [d for d in candidate if d not in fixed_set]
        if newly:
            fixed_set.update(newly)
            added.extend(newly)
            nodes.append(nid)

    return StabilizationReport(
        kind="truss_only_rotations",
        added_dofs=sorted(added),
        node_ids=sorted(set(nodes)),
        dof_names=["RX", "RY", "RZ"],
    )

@dataclass
class LoadedFrame:
    """A frame with loads applied, ready for analysis."""
    frame: Frame
    loads: FrameLoadCase

    settings: FrameAnalysisSettings = field(default_factory=FrameAnalysisSettings)
    
    _nodal_displacements: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _reactions: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(init=False, default_factory=dict)

    analysis_result: Optional[FrameAnalysisResult] = field(init=False, default=None)

    # Auto-splitting / bundling metadata
    original_frame: Frame = field(init=False)
    original_loads: FrameLoadCase = field(init=False)
    _member_bundle: Dict[str, List[str]] = field(init=False, default_factory=dict)  # parent -> ordered segment ids
    _member_segment_parent: Dict[str, str] = field(init=False, default_factory=dict)  # segment -> parent
    _member_segment_offset: Dict[str, float] = field(init=False, default_factory=dict)  # segment -> x0 along parent
    _member_nodes_along: Dict[str, List[Tuple[float, str]]] = field(init=False, default_factory=dict)  # parent -> [(x,node_id)]

    # Diagnostics for batch extraction
    _to_loaded_beams_failures: Dict[str, str] = field(init=False, default_factory=dict)
    
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

        expanded_frame = Frame.from_nodes_and_members(list(nodes_by_id.values()), new_members)

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

        if self.settings.imperfection_model == "notional_loads":
            F_global = F_global + self._build_notional_load_vector(node_to_idx, 6 * len(node_ids))

        fixed_dofs_physical, constrained_nodes = _collect_fixed_dofs(self.frame, node_to_idx)
        stabilization = _truss_only_rotation_stabilization(self.frame, node_to_idx, fixed_dofs_physical)
        fixed_dofs_total = sorted(set(fixed_dofs_physical) | set(stabilization.added_dofs))

        cable_member_ids = [m.id for m in self.frame.members if m.element_type == "cable"]
        member_axial_scales: Dict[str, float] = {}
        for mid in cable_member_ids:
            member_axial_scales[mid] = 1.0

        member_stiffness_scales: Dict[str, ElementStiffnessScales] = {}
        for member in self.frame.members:
            if member.element_type != "beam":
                continue
            member_stiffness_scales[member.id] = ElementStiffnessScales(
                A=self.settings.axial_stiffness_factor,
                Iy=self.settings.bending_stiffness_factor,
                Iz=self.settings.bending_stiffness_factor,
                J=self.settings.torsion_stiffness_factor,
            )

        # Iteration settings
        slack_scale = 1e-9
        tol_n = 1e-6
        max_iter = int(self.settings.max_iter)
        n_steps = int(self.settings.n_steps)
        omega = float(self.settings.relaxation_omega)

        d_global: Optional[np.ndarray] = None
        K_global: Optional[np.ndarray] = None

        warnings: List[str] = []
        design_grade_notes: List[str] = []

        if self.settings.imperfection_model == "notional_loads":
            design_grade_notes.append(
                f"Notional loads enabled (factor={self.settings.notional_factor}, p_source={self.settings.notional_p_source})"
            )

        if (
            abs(float(self.settings.bending_stiffness_factor) - 1.0) > 1e-12
            or abs(float(self.settings.axial_stiffness_factor) - 1.0) > 1e-12
            or abs(float(self.settings.torsion_stiffness_factor) - 1.0) > 1e-12
        ):
            design_grade_notes.append(
                "Element stiffness scaling enabled "
                f"(A={self.settings.axial_stiffness_factor}, Iy/Iz={self.settings.bending_stiffness_factor}, J={self.settings.torsion_stiffness_factor})"
            )

        if self.settings.second_order:
            design_grade_notes.append(
                f"Second-order analysis enabled (method={self.settings.analysis_method}, n_steps={n_steps}, omega={omega})"
            )

        def _assemble_global_stiffness() -> tuple[np.ndarray, dict]:
            if cable_member_ids:
                return assemble_global_stiffness(
                    self.frame,
                    node_to_idx,
                    member_axial_scales=member_axial_scales,
                    member_stiffness_scales=member_stiffness_scales,
                )
            return assemble_global_stiffness(
                self.frame,
                node_to_idx,
                member_stiffness_scales=member_stiffness_scales,
            )

        def _solve_tangent(K: np.ndarray, Kg: np.ndarray, F: np.ndarray, fixed_dofs: List[int]) -> np.ndarray:
            return solve_displacements(K + Kg, F, fixed_dofs)

        def _second_order_solve(
            F_step: np.ndarray,
            member_end_forces_seed: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]],
            d_seed: Optional[np.ndarray],
        ) -> tuple[np.ndarray, dict, Dict[str, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]], bool, int]:
            """Inner second-order solve (fixed cable state).

            Returns:
                (d_global, m_mats, nodal_displacements, member_end_forces, converged, iters)
            """

            last_d = d_seed
            last_member_end_forces = member_end_forces_seed

            for it in range(max_iter):
                K_lin, mats = _assemble_global_stiffness()

                if last_member_end_forces is None:
                    Kg = np.zeros_like(K_lin)
                else:
                    Kg = assemble_geometric_stiffness(
                        self.frame,
                        node_to_idx,
                        last_member_end_forces,
                        member_matrices=mats,
                    )

                d_new = _solve_tangent(K_lin, Kg, F_step, fixed_dofs_total)

                if last_d is not None and abs(omega - 1.0) > 1e-12:
                    d_new = (1.0 - omega) * last_d + omega * d_new

                nodal_displacements, member_end_forces = recover_member_end_forces(
                    self.frame,
                    node_to_idx,
                    d_new,
                    mats,
                )

                if last_d is not None:
                    du = float(np.linalg.norm(d_new - last_d))
                    un = float(np.linalg.norm(d_new))
                    if du <= max(float(self.settings.tol_u_abs), float(self.settings.tol_u_rel) * max(un, 1e-12)):
                        return d_new, mats, nodal_displacements, member_end_forces, True, it + 1

                last_d = d_new
                last_member_end_forces = member_end_forces

            return last_d if last_d is not None else np.zeros(6 * len(node_ids)), mats, nodal_displacements, member_end_forces, False, max_iter

        converged = True
        iterations = 0

        last_member_end_forces: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None
        last_d: Optional[np.ndarray] = None

        # Load stepping (continuation) for second-order robustness.
        for step in range(1, max(n_steps, 1) + 1):
            lam = float(step) / float(max(n_steps, 1))
            F_step = F_global * lam

            # Outer loop: cable tension-only active set (if enabled)
            for cable_it in range(max_iter):
                if self.settings.second_order:
                    d_step, mats, nodal_displacements, member_end_forces, ok, iters = _second_order_solve(
                        F_step,
                        member_end_forces_seed=last_member_end_forces,
                        d_seed=last_d,
                    )
                    iterations += int(iters)
                    if not ok:
                        converged = False
                        warnings.append(
                            f"Second-order iteration did not converge at load step {step}/{n_steps} (cable_outer={cable_it + 1})"
                        )
                else:
                    K_lin, mats = _assemble_global_stiffness()
                    d_step = solve_displacements(K_lin, F_step, fixed_dofs_total)
                    nodal_displacements, member_end_forces = recover_member_end_forces(
                        self.frame,
                        node_to_idx,
                        d_step,
                        mats,
                    )
                    iterations += 1

                # Save for next iteration/step
                d_global = d_step
                last_d = d_step
                last_member_end_forces = member_end_forces

                if not self.settings.cable_tension_only or not cable_member_ids:
                    break

                # Update cable activity
                changed = False
                for mid in cable_member_ids:
                    start_f, _end_f = member_end_forces[mid]
                    n_axial = -float(start_f[0])
                    desired = 1.0 if n_axial >= -tol_n else slack_scale
                    if float(member_axial_scales[mid]) != float(desired):
                        member_axial_scales[mid] = float(desired)
                        changed = True

                if not changed:
                    break

                if cable_it == max_iter - 1:
                    converged = False
                    warnings.append(f"Cable tension-only iteration hit max_iter at load step {step}/{n_steps}")

            # Update final state at this load step
            self._nodal_displacements = nodal_displacements
            self._member_end_forces = member_end_forces

        # Record that stepping/damping were used.
        if self.settings.second_order and n_steps > 1:
            design_grade_notes.append(f"Load stepping enabled (n_steps={n_steps})")
        if self.settings.second_order and abs(omega - 1.0) > 1e-12:
            design_grade_notes.append(f"Under-relaxation enabled (omega={omega})")

        if d_global is None:
            raise RuntimeError("Frame analysis did not run")

        # Recompute reactions from the final linear stiffness (material stiffness only).
        # For second-order, this is an approximation; reaction reporting is primarily for equilibrium checks.
        if K_global is None:
            K_global, _ = _assemble_global_stiffness()

        R_global = K_global @ d_global - F_global
        for nid in constrained_nodes:
            self._reactions[nid] = R_global[node_to_idx[nid]*6 : node_to_idx[nid]*6 + 6]

        if stabilization.used:
            design_grade_notes.append(
                "Numerical stabilization constraints were added (truss/cable-only node rotations)"
            )

        design_grade_ok = bool(converged) and (not stabilization.used)

        self.analysis_result = FrameAnalysisResult(
            nodal_displacements=dict(self._nodal_displacements),
            reactions=dict(self._reactions),
            member_end_forces=dict(self._member_end_forces),
            settings=self.settings,
            converged=bool(converged),
            iterations=int(iterations),
            design_grade_ok=bool(design_grade_ok),
            design_grade_notes=design_grade_notes,
            warnings=warnings,
            stabilization=stabilization,
            fixed_dofs_physical=fixed_dofs_physical,
            fixed_dofs_total=fixed_dofs_total,
        )

    def _build_notional_load_vector(self, node_to_idx: Dict[str, int], n_dofs: int) -> np.ndarray:
        """Build the global notional load vector (imperfection equivalent).

        Default behavior is deterministic and input-based: it uses applied vertical loads
        (nodal forces + member loads lumped to nodes) to compute a gravity measure P, then
        applies lateral notional loads H = alpha * P.
        """

        if self.settings.notional_p_source != "input_based":
            # Reaction-based notional loads are intentionally not implemented yet.
            # It must be explicit because it depends on supports/springs.
            raise NotImplementedError("Reaction-based notional P source is not implemented")

        P_by_node = self._compute_input_based_gravity_p_by_node(node_to_idx)

        F = np.zeros(n_dofs)
        alpha = float(self.settings.notional_factor)
        sx, sy = float(self.settings.notional_signs[0]), float(self.settings.notional_signs[1])

        for nid, P in P_by_node.items():
            H = alpha * float(P)
            base = node_to_idx[nid] * 6
            if "x" in self.settings.notional_axes:
                F[base + 0] += sx * H
            if "y" in self.settings.notional_axes:
                F[base + 1] += sy * H

        return F

    def _compute_input_based_gravity_p_by_node(self, node_to_idx: Dict[str, int]) -> Dict[str, float]:
        """Compute a deterministic gravity-load measure P per node from applied loads.

        Convention:
            - Global Z is treated as "vertical".
            - Downward vertical load is negative Fz.
            - P is returned as a positive scalar, summing only downward contributions.

        This intentionally avoids reactions so P is not sensitive to support stiffness.
        """

        def _downward_p_from_global_force(f_g: np.ndarray) -> float:
            fz = float(f_g[2])
            return max(0.0, -fz)

        P: Dict[str, float] = {nid: 0.0 for nid in self.frame.nodes.keys()}

        # 1) Direct nodal forces
        for nf in self.loads.nodal_forces:
            f_g = nf.force
            if getattr(nf, "coords", "global") == "local":
                m_ref = self._get_reference_member(getattr(nf, "reference_member_id", None))
                f_g = m_ref.transformation_matrix.T @ nf.force
            P[nf.node_id] += _downward_p_from_global_force(f_g)

        # 2) Member point forces (distributed to end nodes like the solver load vector)
        for mpf in self.loads.member_point_forces:
            m = self.frame.get_member(mpf.member_id)
            x = mpf.position * m.length if mpf.position_type == "relative" else mpf.position
            f_g = m.transformation_matrix.T @ mpf.force if mpf.coords == "local" else mpf.force
            p_here = _downward_p_from_global_force(f_g)
            if p_here <= 0.0:
                continue
            N1, N2 = 1 - x / m.length, x / m.length
            P[m.start_node_id] += N1 * p_here
            P[m.end_node_id] += N2 * p_here

        # 3) Member distributed forces (lump total vertical to end nodes)
        for mdf in self.loads.member_distributed_forces:
            m = self.frame.get_member(mdf.member_id)

            # Convert to global if needed
            s_f_g = m.transformation_matrix.T @ mdf.start_force if mdf.coords == "local" else mdf.start_force
            e_f_g = m.transformation_matrix.T @ mdf.end_force if mdf.coords == "local" else mdf.end_force

            # Total resultant (linearly varying): average * length
            Lseg = float(mdf.end_position - mdf.start_position)
            if Lseg <= 0.0:
                continue
            w_avg = 0.5 * (s_f_g + e_f_g)
            F_tot = w_avg * Lseg
            p_here = _downward_p_from_global_force(F_tot)
            if p_here <= 0.0:
                continue

            # Simple deterministic lump: half to each end
            P[m.start_node_id] += 0.5 * p_here
            P[m.end_node_id] += 0.5 * p_here

        # Prune zero entries (but keep determinism by returning only nodes with P>0)
        return {nid: float(v) for nid, v in P.items() if float(v) > 0.0}

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

    @staticmethod
    def _beam_shape_functions(x: float, L: float) -> tuple[float, float, float, float]:
        xi = x / L
        n1 = 1.0 - 3.0 * xi ** 2 + 2.0 * xi ** 3
        n2 = L * (xi - 2.0 * xi ** 2 + xi ** 3)
        n3 = 3.0 * xi ** 2 - 2.0 * xi ** 3
        n4 = L * (-xi ** 2 + xi ** 3)
        return n1, n2, n3, n4

    def _consistent_member_distributed_load(
        self,
        member: Member,
        start_pos: float,
        end_pos: float,
        start_force_l: np.ndarray,
        end_force_l: np.ndarray,
        n_gauss: int = 5,
    ) -> np.ndarray:
        L = member.length
        f_eq = np.zeros(12)
        if end_pos <= start_pos:
            return f_eq

        xi, wi = np.polynomial.legendre.leggauss(n_gauss)
        a = float(start_pos)
        b = float(end_pos)
        half = 0.5 * (b - a)
        mid = 0.5 * (b + a)

        for xg, wg in zip(xi, wi):
            x = mid + half * xg
            if b == a:
                t = 0.0
            else:
                t = (x - a) / (b - a)
            w_l = (1.0 - t) * start_force_l + t * end_force_l

            n1_ax = 1.0 - x / L
            n2_ax = x / L
            n1, n2, n3, n4 = self._beam_shape_functions(x, L)

            weight = wg * half

            if abs(w_l[0]) > 0.0:
                f_eq[0] += n1_ax * w_l[0] * weight
                f_eq[6] += n2_ax * w_l[0] * weight

            if abs(w_l[1]) > 0.0:
                f_eq[1] += n1 * w_l[1] * weight
                f_eq[5] += n2 * w_l[1] * weight
                f_eq[7] += n3 * w_l[1] * weight
                f_eq[11] += n4 * w_l[1] * weight

            if abs(w_l[2]) > 0.0:
                f_eq[2] += n1 * w_l[2] * weight
                f_eq[4] += -n2 * w_l[2] * weight
                f_eq[8] += n3 * w_l[2] * weight
                f_eq[10] += -n4 * w_l[2] * weight

        return f_eq

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
        if mdf.coords == "global":
            s_f_l = m.transformation_matrix @ mdf.start_force
            e_f_l = m.transformation_matrix @ mdf.end_force
        else:
            s_f_l = mdf.start_force
            e_f_l = mdf.end_force

        f_eq_local = self._consistent_member_distributed_load(
            member=m,
            start_pos=float(mdf.start_position),
            end_pos=float(mdf.end_position),
            start_force_l=s_f_l,
            end_force_l=e_f_l,
        )
        T12 = build_transformation_matrix_12x12(m.transformation_matrix)
        f_eq_global = T12.T @ f_eq_local

        s_idx = node_to_idx[m.start_node_id] * 6
        e_idx = node_to_idx[m.end_node_id] * 6
        F[s_idx:s_idx + 6] += f_eq_global[0:6]
        F[e_idx:e_idx + 6] += f_eq_global[6:12]

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

        def _has_any_support(s: Optional[str]) -> bool:
            return bool(s) and any(c == "1" for c in s)

        def _merge_support_types(a: Optional[str], b: Optional[str]) -> Optional[str]:
            if not a:
                return b
            if not b:
                return a
            return "".join("1" if (ca == "1" or cb == "1") else "0" for ca, cb in zip(a, b))

        def _stabilize_supports(supports: list[Support]) -> list[Support]:
            """Ensure the 1D beam FEM has enough boundary conditions.

            The transverse (Hermite) beam solve uses translation+rotation DOFs.
            If *no* Ry/Rz constraints exist anywhere, the stiffness matrix becomes singular.

            Frame models can legitimately omit some rotations/translations (e.g. lifted case),
            so we add a minimal stabilizer at x=0 (or merge into existing x=0 support)
            for any missing DOFs.
            """
            have = [any(s.type[i] == "1" for s in supports) for i in range(6)]

            if all(have):
                return supports

            stabilizer = list("000000")
            for i in range(6):
                if not have[i]:
                    stabilizer[i] = "1"
            stab_type = "".join(stabilizer)

            # Merge into an existing x=0 support if present.
            for s in supports:
                if abs(float(s.x) - 0.0) <= 1e-9:
                    s.type = _merge_support_types(s.type, stab_type) or s.type
                    return supports

            return [Support(x=0.0, type=stab_type)] + supports

        def _supports_from_chain(chain: list[tuple[float, str]]) -> tuple[list[Support], set[str]]:
            """Build 1D supports from frame node constraints along a member chain."""
            # Map x -> merged support type
            x_to_type: dict[float, str] = {}
            supported_node_ids: set[str] = set()

            for x, nid in chain:
                node = self.frame.get_node(nid)
                if not _has_any_support(node.support):
                    continue
                supported_node_ids.add(nid)
                x_f = float(x)
                # Dedup by tolerance
                merged_key = None
                for existing_x in x_to_type.keys():
                    if abs(existing_x - x_f) <= 1e-9:
                        merged_key = existing_x
                        break
                if merged_key is None:
                    x_to_type[x_f] = node.support  # type: ignore[assignment]
                else:
                    x_to_type[merged_key] = _merge_support_types(x_to_type[merged_key], node.support) or x_to_type[merged_key]

            supports = [Support(x=x, type=t) for x, t in sorted(x_to_type.items(), key=lambda kv: kv[0])]
            supports = _stabilize_supports(supports)
            return supports, supported_node_ids

        def _as_global(vec: np.ndarray, coords: str, reference_member_id: Optional[str]) -> np.ndarray:
            if coords != "local":
                return vec
            m_ref = self._get_reference_member(reference_member_id)
            return m_ref.transformation_matrix.T @ vec

        # Reset diagnostics
        self._to_loaded_beams_failures = {}

        # If we didn't split, keep the legacy behavior.
        has_bundles = bool(self._member_bundle) and any(len(v) > 1 for v in self._member_bundle.values())
        if not has_bundles:
            loaded_beams: Dict[str, LoadedBeam] = {}
            for member in self.frame.members:
                try:
                    start_f, end_f = self._member_end_forces[member.id]

                    # Supports from frame node constraints (start/end only in non-bundled path)
                    chain = [(0.0, member.start_node_id), (member.length, member.end_node_id)]
                    supports, supported_node_ids = _supports_from_chain(chain)
                    beam = Beam1D(L=member.length, material=member.material, section=member.section, supports=supports)
                    lc = LoadCase(name=f"Member {member.id} loads")

                    start_supported = member.start_node_id in supported_node_ids
                    end_supported = member.end_node_id in supported_node_ids

                    if member.element_type != "beam":
                        if not end_supported:
                            lc.add_point_force(PointForce(point=np.array([member.length, 0, 0]), force=np.array([-end_f[0], 0, 0])))
                    else:
                        if not end_supported:
                            lc.add_point_force(PointForce(point=np.array([member.length, 0, 0]), force=-end_f[0:3]))
                            if any(abs(end_f[3:6]) > 1e-10):
                                lc.add_moment(Moment(x=member.length, moment=-end_f[3:6]))
                        if not start_supported:
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
                except Exception as e:
                    self._to_loaded_beams_failures[member.id] = str(e)
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

            try:
                chain = self._member_nodes_along.get(parent.id, [])
                supports, supported_node_ids = _supports_from_chain(chain)
                beam = Beam1D(L=parent.length, material=parent.material, section=parent.section, supports=supports)
                lc = LoadCase(name=f"Member {parent.id} loads")
                # Determine whether member ends are supported (by node constraints)
                start_nid = chain[0][1] if chain else parent.start_node_id
                end_nid = chain[-1][1] if chain else parent.end_node_id
                start_supported = start_nid in supported_node_ids
                end_supported = end_nid in supported_node_ids

                if parent.element_type != "beam":
                    if not end_supported:
                        lc.add_point_force(
                            PointForce(
                                point=np.array([parent.length, 0, 0]),
                                force=np.array([-end_f[0], 0, 0]),
                            )
                        )
                else:
                    if not end_supported:
                        lc.add_point_force(PointForce(point=np.array([parent.length, 0, 0]), force=-end_f[0:3]))
                        if any(abs(end_f[3:6]) > 1e-10):
                            lc.add_moment(Moment(x=parent.length, moment=-end_f[3:6]))
                    if not start_supported:
                        if any(abs(start_f[3:6]) > 1e-10):
                            lc.add_moment(Moment(x=0.0, moment=start_f[3:6]))

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
                    lc.add_distributed_force(
                        DistributedForce(
                            np.array([s_x, 0, 0]),
                            np.array([e_x, 0, 0]),
                            s_f_l,
                            e_f_l,
                        )
                    )

                loaded_beams[parent.id] = LoadedBeam(beam, lc)
            except Exception as e:
                self._to_loaded_beams_failures[parent.id] = str(e)

        return loaded_beams

    @property
    def to_loaded_beams_failures(self) -> Dict[str, str]:
        """Per-member errors from the most recent call to to_loaded_beams()."""
        return dict(self._to_loaded_beams_failures)

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

        def _has_any_support(s: Optional[str]) -> bool:
            return bool(s) and any(c == "1" for c in s)

        def _merge_support_types(a: Optional[str], b: Optional[str]) -> Optional[str]:
            if not a:
                return b
            if not b:
                return a
            return "".join("1" if (ca == "1" or cb == "1") else "0" for ca, cb in zip(a, b))

        def _stabilize_supports(supports: list[Support]) -> list[Support]:
            have = [any(s.type[i] == "1" for s in supports) for i in range(6)]
            if all(have):
                return supports
            stabilizer = list("000000")
            for i in range(6):
                if not have[i]:
                    stabilizer[i] = "1"
            stab_type = "".join(stabilizer)
            for s in supports:
                if abs(float(s.x) - 0.0) <= 1e-9:
                    s.type = _merge_support_types(s.type, stab_type) or s.type
                    return supports
            return [Support(x=0.0, type=stab_type)] + supports

        def _supports_from_chain(chain: list[tuple[float, str]]) -> tuple[list[Support], set[str]]:
            x_to_type: dict[float, str] = {}
            supported_node_ids: set[str] = set()
            for x, nid in chain:
                node = self.frame.get_node(nid)
                if not _has_any_support(node.support):
                    continue
                supported_node_ids.add(nid)
                x_f = float(x)
                merged_key = None
                for existing_x in x_to_type.keys():
                    if abs(existing_x - x_f) <= 1e-9:
                        merged_key = existing_x
                        break
                if merged_key is None:
                    x_to_type[x_f] = node.support  # type: ignore[assignment]
                else:
                    x_to_type[merged_key] = _merge_support_types(x_to_type[merged_key], node.support) or x_to_type[merged_key]
            supports = [Support(x=x, type=t) for x, t in sorted(x_to_type.items(), key=lambda kv: kv[0])]
            supports = _stabilize_supports(supports)
            return supports, supported_node_ids
        
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
        
        # Get the node chain for this member and derive supports from frame node constraints.
        chain = self._member_nodes_along.get(member_id, [])
        node_to_x = {nid: float(x) for x, nid in chain}
        supports, supported_node_ids = _supports_from_chain(chain)

        beam = Beam1D(L=parent.length, material=parent.material, section=parent.section, supports=supports)
        lc = LoadCase(name=f"Member {member_id} extracted")

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
        
        # Apply boundary loads at the member end only if that end is NOT modeled as a support.
        # (Otherwise we'd be double-counting the restraint + the end forces.)
        end_nid = parent.end_node_id

        end_supported = (chain[-1][1] if chain else parent.end_node_id) in supported_node_ids

        if (not end_supported) and (end_nid not in fixed_node_ids):
            # Apply end forces and moments as negative (equilibrium)
            lc.add_point_force(PointForce(point=np.array([parent.length, 0, 0]), force=-end_f[0:3]))
            if any(abs(end_f[3:6]) > 1e-10):
                lc.add_moment(Moment(x=parent.length, moment=-end_f[3:6]))
        
        # Check intermediate nodes for reactions and apply as loads only when NOT modeled as supports.
        for x, nid in chain:
            if x <= 1e-9 or x >= parent.length - 1e-9:
                continue  # Skip start/end nodes already handled
            if nid in supported_node_ids:
                continue
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
