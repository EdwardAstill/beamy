from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional, List, Set
import numpy as np

from .frame import Frame
from ..core.loads import (
    FrameLoadCase, MemberPointForce, MemberDistributedForce, 
    PointForce, Moment, DistributedForce, LoadCase
)
from .member import Member
from .solver import analyze_frame_geometry, solve_displacements
from .results import MemberResults
from ..core.math import build_local_stiffness_matrix, build_transformation_matrix_12x12


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

@dataclass
class LoadedFrame:
    """A frame with loads applied, ready for analysis."""
    frame: Frame
    loads: FrameLoadCase
    
    _nodal_displacements: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _reactions: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(init=False, default_factory=dict)
    
    def __post_init__(self) -> None:
        self._validate_loads()
        self._adjust_full_length_loads()
        self._analyze_frame()
    
    def _validate_loads(self) -> None:
        for nf in self.loads.nodal_forces:
            if nf.node_id not in self.frame.nodes: raise ValueError(f"Unknown node '{nf.node_id}'")
        for nm in self.loads.nodal_moments:
            if nm.node_id not in self.frame.nodes: raise ValueError(f"Unknown node '{nm.node_id}'")
        mids = [m.id for m in self.frame.members]
        for mpf in self.loads.member_point_forces:
            if mpf.member_id not in mids: raise ValueError(f"Unknown member '{mpf.member_id}'")
        for mdf in self.loads.member_distributed_forces:
            if mdf.member_id not in mids: raise ValueError(f"Unknown member '{mdf.member_id}'")
    
    def _adjust_full_length_loads(self) -> None:
        for mdf in self.loads.member_distributed_forces:
            if mdf.end_position == -1.0:
                mdf.end_position = self.frame.get_member(mdf.member_id).length
    
    def _analyze_frame(self) -> None:
        node_ids = sorted(self.frame.nodes.keys())
        node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
        K_global, m_mats = analyze_frame_geometry(self.frame, node_to_idx)
        F_global = self._build_load_vector(node_to_idx, 6 * len(node_ids))
        
        fixed_dofs, constrained_nodes = _collect_fixed_dofs(self.frame, node_to_idx)
        
        d_global = solve_displacements(K_global, F_global, fixed_dofs)
        for nid in node_ids: self._nodal_displacements[nid] = d_global[node_to_idx[nid]*6 : node_to_idx[nid]*6 + 6]
        
        R_global = K_global @ d_global - F_global
        for nid in constrained_nodes:
            self._reactions[nid] = R_global[node_to_idx[nid]*6 : node_to_idx[nid]*6 + 6]
            
        for member in self.frame.members:
            d_m_g = np.concatenate([self._nodal_displacements[member.start_node_id], self._nodal_displacements[member.end_node_id]])
            k_l, T = m_mats[member.id]
            f_m_l = k_l @ (T @ d_m_g)
            self._member_end_forces[member.id] = (f_m_l[0:6], f_m_l[6:12])

    def _build_load_vector(self, node_to_idx: Dict[str, int], n_dofs: int) -> np.ndarray:
        F = np.zeros(n_dofs)
        for nf in self.loads.nodal_forces: F[node_to_idx[nf.node_id]*6 : node_to_idx[nf.node_id]*6+3] += nf.force
        for nm in self.loads.nodal_moments: F[node_to_idx[nm.node_id]*6+3 : node_to_idx[nm.node_id]*6+6] += nm.moment
        for mpf in self.loads.member_point_forces: self._add_member_point_force_to_vector(mpf, node_to_idx, F)
        for mdf in self.loads.member_distributed_forces: self._add_member_distributed_force_to_vector(mdf, node_to_idx, F)
        return F

    def _add_member_point_force_to_vector(self, mpf, node_to_idx, F):
        m = self.frame.get_member(mpf.member_id)
        x = mpf.position * m.length if mpf.position_type == "relative" else mpf.position
        f_g = m.transformation_matrix.T @ mpf.force if mpf.coords == "local" else mpf.force
        N1, N2 = 1 - x/m.length, x/m.length
        F[node_to_idx[m.start_node_id]*6 : node_to_idx[m.start_node_id]*6+3] += N1 * f_g
        F[node_to_idx[m.end_node_id]*6 : node_to_idx[m.end_node_id]*6+3] += N2 * f_g

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
        
        loaded_beams = {}
        for member in self.frame.members:
            start_f, end_f = self._member_end_forces[member.id]
            beam = Beam1D(L=member.length, material=member.material, section=member.section, supports=[Support(x=0.0, type="111111")])
            lc = LoadCase(name=f"Member {member.id} loads")
            lc.add_point_force(PointForce(point=np.array([member.length, 0, 0]), force=-end_f[0:3]))
            if any(abs(end_f[3:6]) > 1e-10): lc.add_moment(Moment(x=member.length, moment=-end_f[3:6]))
            for mpf in self.loads.member_point_forces:
                if mpf.member_id == member.id:
                    x = mpf.position * member.length if mpf.position_type == "relative" else mpf.position
                    f_l = member.transformation_matrix @ mpf.force if mpf.coords == "global" else mpf.force
                    lc.add_point_force(PointForce(point=np.array([x, 0, 0]), force=f_l))
            for mdf in self.loads.member_distributed_forces:
                if mdf.member_id == member.id:
                    s_f_l = member.transformation_matrix @ mdf.start_force if mdf.coords == "global" else mdf.start_force
                    e_f_l = member.transformation_matrix @ mdf.end_force if mdf.coords == "global" else mdf.end_force
                    lc.add_distributed_force(DistributedForce(np.array([mdf.start_position, 0, 0]), np.array([mdf.end_position, 0, 0]), s_f_l, e_f_l))
            loaded_beams[member.id] = LoadedBeam(beam, lc)
        return loaded_beams

    def get_member_results(self, member_id: str) -> MemberResults:
        return MemberResults(member=self.frame.get_member(member_id), loaded_beam=self.to_loaded_beams()[member_id])

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
