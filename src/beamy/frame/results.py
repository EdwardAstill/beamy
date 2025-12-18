from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING, Tuple, List
import numpy as np

from ..core.results import Result, AnalysisResult

if TYPE_CHECKING:
    from .member import Member
    from ..beam1d.analysis import LoadedBeam


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
    member_loads: List[Tuple[float, float, np.ndarray]],  # [(start_x, end_x, force_per_m), ...]
    points: int = 100,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute internal forces directly from member end forces using equilibrium.
    
    This avoids the broken 1D FEM cantilever approach which gives wrong reactions.
    
    Args:
        member: The Member object
        start_f: Local forces at start [Fx, Fy, Fz, Mx, My, Mz]
        end_f: Local forces at end [Fx, Fy, Fz, Mx, My, Mz]
        member_loads: List of distributed loads on this member
        points: Number of evaluation points
        
    Returns:
        (xs, N, Vy, Vz, T, My, Mz) - position and internal forces
        
    Sign convention for internal forces:
        - N > 0: tension
        - V > 0: shear in positive local axis direction  
        - M > 0: moment causing positive curvature
    """
    L = member.length
    xs = np.linspace(0, L, points)
    
    # Initialize with forces from start node equilibrium
    # Internal force at cut x = -(sum of external forces from 0 to x)
    # At x=0+: internal = -start_f (member reacts against node force)
    N = np.full(points, -start_f[0])
    Vy = np.full(points, -start_f[1])
    Vz = np.full(points, -start_f[2])
    T = np.full(points, -start_f[3])
    
    # Bending moments vary linearly with distance due to shear
    # My(x) = My(0) + integral of Vz from 0 to x
    # Mz(x) = Mz(0) - integral of Vy from 0 to x (sign from cross product)
    My = -start_f[4] + (-start_f[2]) * xs  # My from Vz acting over arm x
    Mz = -start_f[5] - (-start_f[1]) * xs  # Mz from Vy acting over arm x (negative due to sign convention)
    
    # Add effects of distributed member loads
    for load_start, load_end, w in member_loads:
        # w is force per unit length [wx, wy, wz] in local coords
        for i, x in enumerate(xs):
            if x <= load_start:
                continue
            
            # Length of load acting on segment [load_start, x]
            x_eff = min(x, load_end) - load_start
            if x_eff <= 0:
                continue
                
            # Resultant force from distributed load
            F_dist = w * x_eff
            
            # Position of resultant (midpoint of loaded segment)
            x_res = load_start + x_eff / 2
            arm = x - x_res
            
            # Add to internal forces (internal = -external)
            N[i] -= F_dist[0]
            Vy[i] -= F_dist[1]
            Vz[i] -= F_dist[2]
            My[i] -= F_dist[2] * arm
            Mz[i] += F_dist[1] * arm  # sign convention
    
    return xs, N, Vy, Vz, T, My, Mz


@dataclass
class MemberResultsDirect:
    """
    Direct computation of member internal forces from frame analysis end forces.
    
    This replaces the broken LoadedBeam approach which used a 1D FEM cantilever
    model that gave wrong boundary conditions.
    """
    member: "Member"
    start_f: np.ndarray
    end_f: np.ndarray
    member_loads: List[Tuple[float, float, np.ndarray]]
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


# Keep old class for backward compatibility with 1D beam analysis
@dataclass
class MemberResults:
    """Detailed analysis results for an individual member (legacy 1D approach)."""
    member: "Member"
    loaded_beam: "LoadedBeam"
    
    @property
    def axial(self): return self.loaded_beam.axial()
    @property
    def shear_y(self): return self.loaded_beam.shear("y")
    @property
    def shear_z(self): return self.loaded_beam.shear("z")
    @property
    def bending_y(self): return self.loaded_beam.bending("y")
    @property
    def bending_z(self): return self.loaded_beam.bending("z")
    @property
    def torsion(self): return self.loaded_beam.torsion()
    @property
    def von_mises(self): return self.loaded_beam.von_mises()
