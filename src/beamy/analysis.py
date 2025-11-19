# analysis.py
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .beam import Beam1D
from .loads import LoadCase


import numpy as np

def solve_axial_forces(
    L: float,
    support_positions: list[float],
    point_forces: list[tuple[float, float]],            # (x, Fx)
    dist_forces: list[tuple[float, float, float]]       # (x_start, x_end, n_x)
):
    """
    Compute axial reactions and internal axial force N(x).

    Args:
        L : float
            Length of the beam.
        support_positions : list[float]
            x-positions of supports that restrain Ux.
        point_forces : list[(x, Fx)]
            Axial point forces (positive = tension).
        dist_forces : list[(x_start, x_end, n_x)]
            Axial distributed loads in N/m (positive = tension).

    Returns:
        reactions : dict {x_support : Rx}
        N(x) : callable internal axial force along beam (positive = tension)
    """

    # ------------------------------------------
    # 1. SORT INPUTS
    # ------------------------------------------
    support_positions = sorted(support_positions)
    point_forces = sorted(point_forces, key=lambda p: p[0])
    dist_forces = sorted(dist_forces, key=lambda d: d[0])

    n_supports = len(support_positions)

    if n_supports == 0:
        raise ValueError("Beam has no axial restraints (no Ux support). Unstable.")
    
    # ------------------------------------------
    # 2. ASSEMBLE TOTAL EXTERNAL AXIAL FORCE
    # ------------------------------------------
    total_F_points = sum(F for _, F in point_forces)
    total_F_dist = sum(n_x * (x_end - x_start) for x_start, x_end, n_x in dist_forces)
    total_F = total_F_points + total_F_dist

    # ------------------------------------------
    # 3. DETERMINE REACTIONS AT SUPPORTS
    # ------------------------------------------
    # Case A: Exactly 1 axial support → statically determinate
    if n_supports == 1:
        x_s = support_positions[0]
        reactions = {x_s: -total_F}

    # Case B: ≥2 supports → statically indeterminate
    # For v1: distribute reactions proportionally by stiffness of "segments"
    # (this is identical to a 1D axial bar stiffness method with identical EA)
    else:
        # All axial elements assumed identical: stiffness proportional to length⁻¹
        # Build element lengths between supports
        lengths = np.diff(support_positions)
        
        # Build stiffness matrix for a chain of springs (bar elements)
        k_list = [1.0 / L_i for L_i in lengths]

        # Global stiffness matrix
        K = np.zeros((n_supports, n_supports))
        for i, k in enumerate(k_list):
            K[i, i]     += k
            K[i, i + 1] += -k
            K[i + 1, i] += -k
            K[i + 1, i + 1] += k

        # Equivalent nodal loads = external loads collapsed to nearest supports
        # v1: simplest possible consistent load: distribute proportionally by distance.
        F_support = np.zeros(n_supports)

        for x_p, Fp in point_forces:
            # Find nearest support on left
            idx = np.searchsorted(support_positions, x_p) - 1
            if idx < 0:
                # Before first support
                F_support[0] += Fp
            elif idx >= n_supports - 1:
                # After last support
                F_support[-1] += Fp
            else:
                # Between support idx and idx+1 → distribute by distance
                xL = support_positions[idx]
                xR = support_positions[idx + 1]
                t = (x_p - xL) / (xR - xL)
                F_support[idx]     += Fp * (1 - t)
                F_support[idx + 1] += Fp * t

        # Distributed forces → equivalent nodal loads
        for x_start, x_end, n_x in dist_forces:
            # Use centroid of the distribution
            x_c = 0.5 * (x_start + x_end)
            F_total = n_x * (x_end - x_start)

            idx = np.searchsorted(support_positions, x_c) - 1
            if idx < 0:
                F_support[0] += F_total
            elif idx >= n_supports - 1:
                F_support[-1] += F_total
            else:
                xL = support_positions[idx]
                xR = support_positions[idx + 1]
                t = (x_c - xL) / (xR - xL)
                F_support[idx]     += F_total * (1 - t)
                F_support[idx + 1] += F_total * t

        # Solve K * R = -F_support
        R = np.linalg.solve(K, -F_support)

        reactions = {x_s: R[i] for i, x_s in enumerate(support_positions)}

    # ------------------------------------------
    # 4. BUILD INTERNAL FORCE FUNCTION N(x)
    # ------------------------------------------
    def N(x: float) -> float:
        """
        Internal axial force at position x (positive = tension).
        """
        if x < 0 or x > L:
            return 0.0

        # Start with reactions on the left side of the cut
        total_left = 0.0

        for x_s, Rs in reactions.items():
            if x_s <= x:
                total_left += Rs

        # Add point forces on the left
        for x_p, Fp in point_forces:
            if x_p <= x:
                total_left += Fp

        # Add distributed forces integrated on the left
        for x_start, x_end, n_x in dist_forces:
            if x <= x_start:
                continue
            x_eff = min(x, x_end) - x_start
            if x_eff > 0:
                total_left += n_x * x_eff

        # Internal force balances external
        return -total_left

    return reactions, N




@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase

    def analyze(self, points: int = 100):

        #firstly calculate axial forces
        #calculate the axial forces at the ends



        #for directions z and y
            # calculate the shear forces

            #then calculate the bending moments (by also looking at the moment components)

            #then calculate the deflections in each direction does this change when it is combined what if assuming linearity???

        # then use teh x component of the moment vectors to calculation the torsion

