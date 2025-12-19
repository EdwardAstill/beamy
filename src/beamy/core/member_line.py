"""
Shared member line demand computation.

This module provides the canonical equilibrium-based reconstruction of internal
actions along a member, given:
  - Local end forces at start/end
  - Distributed loads along the member (linearly varying)

Used by both:
  - Frame analysis (MemberDemandProvider)
  - 1D analysis (LoadedMember wrapping to frame backend)

This ensures consistent results regardless of analysis path.
"""
from __future__ import annotations
from typing import List, Tuple, TYPE_CHECKING
import numpy as np
import warnings

if TYPE_CHECKING:
    pass


def compute_member_actions(
    length: float,
    start_forces: np.ndarray,
    end_forces: np.ndarray,
    distributed_loads: List[Tuple[float, float, np.ndarray, np.ndarray]],
    points: int = 100,
    check_equilibrium: bool = True,
    equilibrium_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute internal forces along a member from end forces using equilibrium.
    
    This is the canonical implementation used by both frame and 1D analysis.
    
    Args:
        length: Member length
        start_forces: Local forces at start [Fx, Fy, Fz, Mx, My, Mz]
        end_forces: Local forces at end [Fx, Fy, Fz, Mx, My, Mz]
        distributed_loads: List of [(start_x, end_x, start_force, end_force), ...]
                          where forces are [wx, wy, wz] per unit length in local coords
        points: Number of evaluation points along member
        check_equilibrium: If True, verify element equilibrium (warn on mismatch)
        equilibrium_tol: Relative tolerance for equilibrium check
        
    Returns:
        (xs, N, Vy, Vz, T, My, Mz) - positions and internal action components
        
    Sign convention for internal forces:
        - N > 0: tension (pulling member apart)
        - V > 0: shear in positive local axis direction
        - T > 0: torsion about positive local x-axis
        - M > 0: moment causing positive curvature (tension on positive face)
    """
    L = float(length)
    xs = np.linspace(0, L, points)
    
    # Initialize with forces from start node equilibrium.
    # Internal force at a cut just to the right of x=0 must balance the external
    # force applied at the start node. By convention, the internal force
    # is the negative of the node force (action-reaction).
    N = np.full(points, -start_forces[0])
    Vy = np.full(points, -start_forces[1])
    Vz = np.full(points, -start_forces[2])
    T = np.full(points, -start_forces[3])
    
    # Bending moments vary with distance due to shear and applied moments.
    # My(x) is driven by Vz (shear in z-direction creates moment about y)
    # Mz(x) is driven by Vy (shear in y-direction creates moment about z)
    # Starting moment at x=0: internal = -external (reaction to applied moment)
    # Then integrate: dM/dx = V (with appropriate signs)
    My = -start_forces[4] + (-start_forces[2]) * xs
    Mz = -start_forces[5] - (-start_forces[1]) * xs
    
    # Track total distributed load for equilibrium check
    F_dist_total = np.zeros(3)
    M_dist_total_about_start = np.zeros(3)
    
    # Add effects of distributed loads
    for load_start, load_end, w_start, w_end in distributed_loads:
        load_len = load_end - load_start
        if load_len <= 0:
            continue
        
        # Total resultant of this distributed load segment
        F_seg = (w_start + w_end) * load_len / 2.0
        F_dist_total += F_seg
        
        # Moment arm: for linearly varying load, centroid is at:
        # x_c = a + L*(w0 + 2*w1) / (3*(w0+w1))  [per component; fallback to midpoint]
        denom_full = (w_start + w_end)
        x_c_full = np.full(3, (load_start + load_end) / 2.0)
        mask_full = np.abs(denom_full) > 1e-12
        x_c_full[mask_full] = (
            load_start
            + load_len * (w_start[mask_full] + 2.0 * w_end[mask_full]) / (3.0 * denom_full[mask_full])
        )
        
        # Moment about start: r x F, with r=[x_c, 0, 0] => M=[0, -x_c*Fz, x_c*Fy]
        M_dist_total_about_start[1] += -x_c_full[2] * F_seg[2]
        M_dist_total_about_start[2] += x_c_full[1] * F_seg[1]
        
        # For each evaluation point, add the effect of the load up to that point
        for i, x in enumerate(xs):
            if x <= load_start:
                continue
            
            # Effective length: how much of this load segment acts before x
            x_eff = min(x, load_end) - load_start
            if x_eff <= 0:
                continue
            
            # Linearly interpolated load intensity at x_eff
            t_eff = x_eff / load_len if load_len > 0 else 0.0
            w_at_xeff = w_start + (w_end - w_start) * t_eff
            
            # Resultant force from load_start to load_start+x_eff (trapezoidal)
            F_dist = (w_start + w_at_xeff) * x_eff / 2.0
            
            # Centroid of this trapezoid (per component)
            denom = (w_start + w_at_xeff)
            x_c = np.full(3, x_eff / 2.0)
            mask = np.abs(denom) > 1e-12
            x_c[mask] = x_eff * (w_start[mask] + 2.0 * w_at_xeff[mask]) / (3.0 * denom[mask])
            
            # Position of the centroid measured from member start
            x_res = load_start + x_c
            
            # Moment arm from centroid to the cut at x
            arm = x - x_res
            
            # Internal forces = -(external loads up to this point)
            N[i] -= F_dist[0]
            Vy[i] -= F_dist[1]
            Vz[i] -= F_dist[2]
            
            # Moments: My from wz, Mz from wy
            # Sign convention: positive moment creates positive curvature
            My[i] -= F_dist[2] * arm[2]  # wz causes My
            Mz[i] += F_dist[1] * arm[1]  # wy causes Mz (note sign)
    
    # Equilibrium check: verify element free-body balance
    if check_equilibrium:
        # Sum of forces: start_F + end_F + F_dist_total should equal zero
        force_residual = start_forces[:3] + end_forces[:3] + F_dist_total
        
        # Sum of moments about start:
        # start_M + end_M + (r_end x end_F) + M_dist_total should equal zero
        r_end = np.array([L, 0.0, 0.0])
        m_from_end_force = np.cross(r_end, end_forces[:3])
        moment_residual = (
            start_forces[3:] + end_forces[3:] + m_from_end_force + M_dist_total_about_start
        )
        
        # Compute relative residual norm
        scale = max(1.0, float(np.linalg.norm(np.concatenate([start_forces, end_forces]))))
        resid_norm = float(np.linalg.norm(np.concatenate([force_residual, moment_residual])))
        
        if resid_norm / scale > equilibrium_tol:
            warnings.warn(
                f"Member equilibrium residual {resid_norm/scale:.2%}. "
                f"This may indicate inconsistent load/end-force data.",
                RuntimeWarning,
            )
    
    return xs, N, Vy, Vz, T, My, Mz
