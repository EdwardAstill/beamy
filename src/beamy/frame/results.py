from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING, Tuple, List, Dict, Optional
import numpy as np

from ..core.results import Result, AnalysisResult
from ..core.member_line import compute_member_actions

if TYPE_CHECKING:
    from .member import Member
    from ..beam1d.analysis import LoadedMember
    from ..core.material import Material
    from sectiony import Section
    from .frame import Frame
    from ..core.loads import FrameLoadCase


class _ActionOnlyResult:
    """Wrapper to provide .action attribute for compatibility with AnalysisResult interface."""
    def __init__(self, result: Result):
        self.action = result
        # Also forward common Result attributes for direct access
        self._result = result
    
    def at(self, x: float) -> float:
        return self._result.at(x)
    
    @property
    def max(self) -> float:
        return self._result.max
    
    @property
    def min(self) -> float:
        return self._result.min
    
    @property
    def abs_max(self) -> float:
        return self._result.abs_max


def _compute_internal_forces_direct(
    member: "Member",
    start_f: np.ndarray,
    end_f: np.ndarray,
    member_loads: List[Tuple[float, float, np.ndarray, np.ndarray]],
    points: int = 100,
    check_equilibrium: bool = True,
    equilibrium_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute internal forces directly from member end forces using equilibrium.
    
    This delegates to the shared member_line module for consistent behavior
    across 1D and 3D analysis paths.
    
    Args:
        member: The Member object
        start_f: Local forces at start [Fx, Fy, Fz, Mx, My, Mz]
        end_f: Local forces at end [Fx, Fy, Fz, Mx, My, Mz]
        member_loads: List of distributed loads [(start_x, end_x, start_force, end_force), ...]
                      where forces are [wx, wy, wz] per unit length in local coords
        points: Number of evaluation points
        check_equilibrium: If True, verify element equilibrium (warn on mismatch)
        equilibrium_tol: Relative tolerance for equilibrium check
        
    Returns:
        (xs, N, Vy, Vz, T, My, Mz) - position and internal forces
        
    Sign convention for internal forces:
        - N > 0: tension
        - V > 0: shear in positive local axis direction  
        - M > 0: moment causing positive curvature
    """
    return compute_member_actions(
        length=member.length,
        start_forces=start_f,
        end_forces=end_f,
        distributed_loads=member_loads,
        points=points,
        check_equilibrium=check_equilibrium,
        equilibrium_tol=equilibrium_tol,
    )


@dataclass
class MemberResultsDirect:
    """
    Direct computation of member internal forces from frame analysis end forces.
    
    This replaces the broken 1D-only approach which used a 1D FEM cantilever
    model that gave wrong boundary conditions.
    
    member_loads format: [(start_x, end_x, start_force, end_force), ...]
    where start_force and end_force are [wx, wy, wz] per unit length in local coords.
    """
    member: "Member"
    start_f: np.ndarray
    end_f: np.ndarray
    member_loads: List[Tuple[float, float, np.ndarray, np.ndarray]]
    _cache: dict = None
    
    def __post_init__(self):
        self._cache = {}
    
    def _compute(self, points: int = 100):
        key = points
        if key not in self._cache:
            xs, N, Vy, Vz, T, My, Mz = _compute_internal_forces_direct(
                self.member, self.start_f, self.end_f, self.member_loads, points
            )
            self._cache[key] = (xs, N, Vy, Vz, T, My, Mz)
        return self._cache[key]
    
    @property
    def axial(self) -> _ActionOnlyResult:
        xs, N, *_ = self._compute()
        return _ActionOnlyResult(Result(xs, N))
    
    @property
    def shear_y(self) -> _ActionOnlyResult:
        xs, _, Vy, *_ = self._compute()
        return _ActionOnlyResult(Result(xs, Vy))
    
    @property
    def shear_z(self) -> _ActionOnlyResult:
        xs, _, _, Vz, *_ = self._compute()
        return _ActionOnlyResult(Result(xs, Vz))
    
    @property
    def torsion(self) -> _ActionOnlyResult:
        xs, _, _, _, T, *_ = self._compute()
        return _ActionOnlyResult(Result(xs, T))
    
    @property
    def bending_y(self) -> _ActionOnlyResult:
        xs, _, _, _, _, My, _ = self._compute()
        return _ActionOnlyResult(Result(xs, My))
    
    @property
    def bending_z(self) -> _ActionOnlyResult:
        xs, _, _, _, _, _, Mz = self._compute()
        return _ActionOnlyResult(Result(xs, Mz))
    
    @property
    def von_mises(self) -> _ActionOnlyResult:
        """Compute von Mises stress from internal forces."""
        xs, N, Vy, Vz, T, My, Mz = self._compute()
        
        # Section properties
        A = self.member.section.A
        Iy = self.member.section.Iy
        Iz = self.member.section.Iz
        J = self.member.section.J
        y_max = self.member.section.y_max
        z_max = self.member.section.z_max
        r_max = max(abs(y_max), abs(z_max))
        
        # Axial stress: sigma_axial = N / A
        sigma_axial = N / A if A > 0 else np.zeros_like(N)
        
        # Bending stress (max at extreme fiber):
        # sigma_bending = My * z_max / Iy + Mz * y_max / Iz
        sigma_bending_y = np.abs(My) * z_max / Iy if Iy > 0 else np.zeros_like(My)
        sigma_bending_z = np.abs(Mz) * y_max / Iz if Iz > 0 else np.zeros_like(Mz)
        
        # Total normal stress (simplified: max of combined)
        sigma = np.abs(sigma_axial) + sigma_bending_y + sigma_bending_z
        
        # Shear stress (approximate: V/A for shear, T*r/J for torsion)
        tau_shear = (np.abs(Vy) + np.abs(Vz)) / A if A > 0 else np.zeros_like(Vy)
        tau_torsion = np.abs(T) * r_max / J if J > 0 else np.zeros_like(T)
        tau = tau_shear + tau_torsion
        
        # Von Mises: sqrt(sigma^2 + 3*tau^2)
        vm = np.sqrt(sigma**2 + 3 * tau**2)
        
        return _ActionOnlyResult(Result(xs, vm))


@dataclass
class MemberResults:
    """Detailed analysis results for an individual member (legacy 1D approach)."""
    member: "Member"
    loaded_member: "LoadedMember"
    
    @property
    def axial(self): return self.loaded_member.axial()
    @property
    def shear_y(self): return self.loaded_member.shear("y")
    @property
    def shear_z(self): return self.loaded_member.shear("z")
    @property
    def bending_y(self): return self.loaded_member.bending("y")
    @property
    def bending_z(self): return self.loaded_member.bending("z")
    @property
    def torsion(self): return self.loaded_member.torsion()
    @property
    def von_mises(self): return self.loaded_member.von_mises()


@dataclass
class MemberActionProfile:
    """Snapshot of member actions along its length (frame-derived, no re-solve)."""

    member_id: str
    length: float
    material: "Material"
    section: "Section"
    axial: Result
    shear_y: Result
    shear_z: Result
    torsion: Result
    bending_y: Result
    bending_z: Result

    def envelopes(self) -> Dict[str, Tuple[float, float]]:
        return {
            "N": (float(self.axial.min), float(self.axial.max)),
            "Vy": (float(self.shear_y.min), float(self.shear_y.max)),
            "Vz": (float(self.shear_z.min), float(self.shear_z.max)),
            "T": (float(self.torsion.min), float(self.torsion.max)),
            "My": (float(self.bending_y.min), float(self.bending_y.max)),
            "Mz": (float(self.bending_z.min), float(self.bending_z.max)),
        }


class MemberDemandProvider:
    """Provides member end forces/actions directly from frame analysis results.
    
    Supports both segment-level and original-member queries when members are split
    due to intermediate nodes or loads.
    """

    def __init__(
        self,
        frame: "Frame",
        loads: "FrameLoadCase",
        member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]],
        member_bundle: Optional[Dict[str, List[str]]] = None,
    ) -> None:
        self._frame = frame
        self._loads = loads
        self._member_end_forces = member_end_forces
        # member_bundle maps original_member_id -> [segment_ids] in order
        self._member_bundle = member_bundle or {}

    def original_member_ids(self) -> List[str]:
        """Return list of original (unsplit) member IDs."""
        if self._member_bundle:
            return list(self._member_bundle.keys())
        return list(self._member_end_forces.keys())

    def segment_ids(self, original_member_id: str) -> List[str]:
        """Return segment IDs for an original member (in order along member)."""
        return self._member_bundle.get(original_member_id, [original_member_id])

    def _distributed_loads_local(self, member: "Member") -> List[Tuple[float, float, np.ndarray, np.ndarray]]:
        """Get distributed loads for a member in local coordinates.
        
        Returns list of (start_x, end_x, start_force, end_force) tuples.
        Supports linearly varying loads.
        """
        member_loads: List[Tuple[float, float, np.ndarray, np.ndarray]] = []
        for mdf in self._loads.member_distributed_forces:
            if mdf.member_id != member.id:
                continue
            if mdf.coords == "global":
                w_start_local = member.transformation_matrix @ mdf.start_force
                w_end_local = member.transformation_matrix @ mdf.end_force
            else:
                w_start_local = mdf.start_force
                w_end_local = mdf.end_force
            member_loads.append((mdf.start_position, mdf.end_position, w_start_local, w_end_local))
        return member_loads

    def end_forces_local(self, member_id: str) -> Tuple[np.ndarray, np.ndarray]:
        """Return local end forces for a segment member.

        If an original (pre-split) member ID is provided, returns the first-segment
        start forces and last-segment end forces.
        """
        if member_id in self._member_end_forces:
            return self._member_end_forces[member_id]

        seg_ids = self._member_bundle.get(member_id)
        if seg_ids:
            start_f, _ = self._member_end_forces[seg_ids[0]]
            _, end_f = self._member_end_forces[seg_ids[-1]]
            return start_f, end_f

        raise KeyError(member_id)

    def end_forces_local_segments(self, original_member_id: str) -> List[Tuple[str, Tuple[np.ndarray, np.ndarray]]]:
        """Return per-segment local end forces for an original member."""
        seg_ids = self._member_bundle.get(original_member_id, [original_member_id])
        return [(seg_id, self._member_end_forces[seg_id]) for seg_id in seg_ids]

    def actions(self, member_id: str, points: int = 201) -> MemberActionProfile:
        # Accept original IDs transparently when bundling is active.
        if member_id not in self._member_end_forces and member_id in self._member_bundle:
            # Convert a total-points request to a reasonable per-segment density.
            seg_ids = self._member_bundle.get(member_id, [member_id])
            points_per_segment = max(11, int(np.ceil(points / max(1, len(seg_ids)))))
            return self.actions_original(member_id, points_per_segment=points_per_segment)

        member = self._frame.get_member(member_id)
        start_f, end_f = self._member_end_forces[member_id]
        mres = MemberResultsDirect(
            member=member,
            start_f=start_f,
            end_f=end_f,
            member_loads=self._distributed_loads_local(member),
        )
        xs, N, Vy, Vz, T, My, Mz = mres._compute(points=points)
        return MemberActionProfile(
            member_id=member_id,
            length=member.length,
            material=member.material,
            section=member.section,
            axial=Result(xs, N),
            shear_y=Result(xs, Vy),
            shear_z=Result(xs, Vz),
            torsion=Result(xs, T),
            bending_y=Result(xs, My),
            bending_z=Result(xs, Mz),
        )

    def envelopes(self, member_id: str, points: int = 201) -> Dict[str, Tuple[float, float]]:
        profile = self.actions(member_id, points=points)
        return profile.envelopes()

    def actions_original(self, original_member_id: str, points_per_segment: int = 101) -> MemberActionProfile:
        """Get action profile for an original (unsplit) member by combining segments.
        
        This stitches together segments that were created when the original member
        was split due to intermediate nodes or point loads.
        
        Args:
            original_member_id: The ID of the original member (before splitting)
            points_per_segment: Number of evaluation points per segment
            
        Returns:
            MemberActionProfile with continuous x-coordinates along original member length
        """
        seg_ids = self._member_bundle.get(original_member_id, [original_member_id])
        
        if len(seg_ids) == 1:
            # No splitting occurred, just return the single segment
            return self.actions(seg_ids[0], points=points_per_segment)
        
        # Collect data from each segment
        all_xs: List[np.ndarray] = []
        all_N: List[np.ndarray] = []
        all_Vy: List[np.ndarray] = []
        all_Vz: List[np.ndarray] = []
        all_T: List[np.ndarray] = []
        all_My: List[np.ndarray] = []
        all_Mz: List[np.ndarray] = []
        
        x_offset = 0.0
        total_length = 0.0
        first_member = None
        
        for i, seg_id in enumerate(seg_ids):
            member = self._frame.get_member(seg_id)
            if first_member is None:
                first_member = member
            
            start_f, end_f = self._member_end_forces[seg_id]
            mres = MemberResultsDirect(
                member=member,
                start_f=start_f,
                end_f=end_f,
                member_loads=self._distributed_loads_local(member),
            )
            xs, N, Vy, Vz, T, My, Mz = mres._compute(points=points_per_segment)
            
            # Offset x-coordinates for continuity
            xs_offset = xs + x_offset
            
            # Avoid duplicating the junction point (except for first segment)
            if i > 0 and len(all_xs) > 0:
                # Skip first point of this segment (it's same location as last point of prev)
                xs_offset = xs_offset[1:]
                N = N[1:]
                Vy = Vy[1:]
                Vz = Vz[1:]
                T = T[1:]
                My = My[1:]
                Mz = Mz[1:]
            
            all_xs.append(xs_offset)
            all_N.append(N)
            all_Vy.append(Vy)
            all_Vz.append(Vz)
            all_T.append(T)
            all_My.append(My)
            all_Mz.append(Mz)
            
            x_offset += member.length
            total_length += member.length
        
        # Concatenate all segments
        xs_combined = np.concatenate(all_xs)
        N_combined = np.concatenate(all_N)
        Vy_combined = np.concatenate(all_Vy)
        Vz_combined = np.concatenate(all_Vz)
        T_combined = np.concatenate(all_T)
        My_combined = np.concatenate(all_My)
        Mz_combined = np.concatenate(all_Mz)
        
        return MemberActionProfile(
            member_id=original_member_id,
            length=total_length,
            material=first_member.material,
            section=first_member.section,  # Assumes uniform section along original member
            axial=Result(xs_combined, N_combined),
            shear_y=Result(xs_combined, Vy_combined),
            shear_z=Result(xs_combined, Vz_combined),
            torsion=Result(xs_combined, T_combined),
            bending_y=Result(xs_combined, My_combined),
            bending_z=Result(xs_combined, Mz_combined),
        )

    def envelopes_original(self, original_member_id: str, points_per_segment: int = 101) -> Dict[str, Tuple[float, float]]:
        """Get force envelopes for an original (unsplit) member."""
        profile = self.actions_original(original_member_id, points_per_segment=points_per_segment)
        return profile.envelopes()
