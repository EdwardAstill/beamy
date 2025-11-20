# analysis.py

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np

from .beam import Beam1D, Support
from .loads import LoadCase


# -------------------------------------------------
# Axial reactions (Fx)
# -------------------------------------------------
def solve_x_forces(
    supports: Sequence[Support],
    loads: LoadCase,
) -> None:
    """
    Compute axial reactions at supports (Fx) for a 1D beam with only axial
    point loads. Modifies support.reactions['Fx'] in place and returns None.

    Args:
        supports:
            Sequence of Support objects. Axial displacement (Ux) is assumed
            to be restrained at these x-positions.
        loads:
            LoadCase. Only the axial component Fx is used here
            (positive = tension along +x).

    Notes:
        - If there is exactly one axial support, the problem is statically
          determinate and the reaction is just minus the sum of axial forces.
        - If there are two or more supports, a 1D axial stiffness chain is
          used to distribute reactions based on relative segment stiffness
          (k ~ 1/L between supports). EA is factored out.

        Assumes loads.convert_dforces_to_pforces() has already been called,
        so all distributed loads are represented as PointForce objects and
        loads.Fxs returns a list of (x, Fx) pairs.
    """
    if not supports:
        raise ValueError("Beam has no axial supports (no Ux restraint). Unstable.")

    # Reset existing axial reactions
    for s in supports:
        s.reactions["Fx"] = 0.0

    # Sort supports by x
    supports_sorted = sorted(supports, key=lambda s: s.x)
    n_supports = len(supports_sorted)

    # List of (x, Fx) from LoadCase
    pf_list: list[tuple[float, float]] = loads.Fxs
    # Total external axial force from point loads
    total_F = sum(Fx for _, Fx in pf_list)

    # -------------------------------------------------
    # Case A: exactly one support → statically determinate
    # -------------------------------------------------
    if n_supports == 1:
        supports_sorted[0].reactions["Fx"] = -total_F
        return

    # -------------------------------------------------
    # Case B: 2+ supports → statically indeterminate
    # use 1D axial stiffness chain with k_i ∝ 1/L_i
    # -------------------------------------------------
    x_supports = np.array([s.x for s in supports_sorted], dtype=float)

    # Segment lengths between supports
    lengths = np.diff(x_supports)
    if np.any(lengths <= 0.0):
        raise ValueError("Support positions must be strictly increasing in x.")

    # Element stiffnesses (EA factored out → k_i = 1/L_i)
    k_list = [1.0 / L_i for L_i in lengths]

    # Global stiffness matrix K (n_supports × n_supports)
    K = np.zeros((n_supports, n_supports), dtype=float)
    for i, k in enumerate(k_list):
        K[i, i]         += k
        K[i, i + 1]     -= k
        K[i + 1, i]     -= k
        K[i + 1, i + 1] += k

    # Equivalent nodal loads at supports (collapse point loads)
    F_support = np.zeros(n_supports, dtype=float)

    for x_p, Fp in pf_list:
        idx = np.searchsorted(x_supports, x_p) - 1

        if idx < 0:
            F_support[0] += Fp
        elif idx >= n_supports - 1:
            F_support[-1] += Fp
        else:
            xL = x_supports[idx]
            xR = x_supports[idx + 1]
            t = (x_p - xL) / (xR - xL)
            F_support[idx]     += Fp * (1.0 - t)
            F_support[idx + 1] += Fp * t

    # Solve K * R = -F_support for reactions R at each support
    R = np.linalg.solve(K, -F_support)

    # Write reactions back into the Support objects (Fx only)
    for i, s in enumerate(supports_sorted):
        s.reactions["Fx"] = float(R[i])


# -------------------------------------------------
# Transverse reactions (Fy/Fz + My/Mz)
# -------------------------------------------------
def solve_transverse_reactions(beam: Beam1D, loads: LoadCase, axis: str = "z") -> None:
    """
    Solve for transverse reactions (shear + bending moment) in one bending plane
    using Euler–Bernoulli 4-DOF beam elements.

    This fills support.reactions[...] *in place* and returns None.

    Args:
        beam : Beam1D
            Beam with .supports (list of Support), .material, .section.
        loads : LoadCase
            Contains applied loads; assumed that dist_forces have already been
            converted to point forces via convert_dforces_to_pforces().
        axis : {"y", "z"}
            Bending plane:
                "z" -> bending in x-z plane (shear Fz, moment My)
                "y" -> bending in x-y plane (shear Fy, moment Mz)

    Notes:
        Uses LoadCase.Fys/Fzs for shear and LoadCase.Mys/Mzs for bending
        (which already include eccentric axial-load moments + explicit moments).
    """
    if axis not in ("y", "z"):
        raise ValueError(f"axis must be 'y' or 'z', got {axis!r}")

    # ---- Map axis to DOF bits / inertia & reaction keys / load lists ----
    if axis == "z":
        # bending in x-z plane, about local y
        trans_bit = 2            # Uz (support.type[2])
        rot_bit = 4              # Ry (support.type[4])
        shear_key = "Fz"
        moment_key = "My"
        I_attr = "Iy"            # beam.section.Iy
        shear_pairs = loads.Fzs  # [(x, Fz)]
        moment_pairs = loads.Mys # [(x, My)]
    else:  # axis == "y"
        # bending in x-y plane, about local z
        trans_bit = 1            # Uy (support.type[1])
        rot_bit = 5              # Rz (support.type[5])
        shear_key = "Fy"
        moment_key = "Mz"
        I_attr = "Iz"            # beam.section.Iz
        shear_pairs = loads.Fys  # [(x, Fy)]
        moment_pairs = loads.Mzs # [(x, Mz)]

    # Get bending stiffness EI
    E = beam.material.E
    I = getattr(beam.section, I_attr, None)
    if I is None:
        raise AttributeError(
            f"Section has no attribute {I_attr!r} required for bending about this axis."
        )
    EI = E * I

    # ---- Collect and sort supports (nodes) ----
    if not beam.supports:
        raise ValueError("Beam has no supports defined.")

    supports_sorted: list[Support] = sorted(beam.supports, key=lambda s: s.x)
    n_nodes = len(supports_sorted)
    dof_per_node = 2
    ndof = dof_per_node * n_nodes

    # Reset reactions in this bending plane
    for s in supports_sorted:
        s.reactions[shear_key] = 0.0
        s.reactions[moment_key] = 0.0

    # Node x-coordinates
    x_nodes = np.array([s.x for s in supports_sorted], dtype=float)

    # ---- Build global stiffness matrix K (ndof x ndof) ----
    K = np.zeros((ndof, ndof), dtype=float)

    # Loop over elements between successive supports
    for e in range(n_nodes - 1):
        x1 = x_nodes[e]
        x2 = x_nodes[e + 1]
        L = x2 - x1
        if L <= 0.0:
            raise ValueError("Support positions must be strictly increasing in x.")

        # Euler–Bernoulli beam element stiffness matrix (local DOF order: w1, θ1, w2, θ2)
        L2 = L * L
        L3 = L2 * L
        k = EI / L3
        ke = k * np.array([
            [ 12.0,     6.0*L,   -12.0,     6.0*L],
            [  6.0*L,   4.0*L2,   -6.0*L,   2.0*L2],
            [-12.0,    -6.0*L,    12.0,    -6.0*L],
            [  6.0*L,   2.0*L2,   -6.0*L,   4.0*L2],
        ])

        # Global DOF indices for this element
        i_w = dof_per_node * e
        i_t = i_w + 1
        j_w = dof_per_node * (e + 1)
        j_t = j_w + 1
        idx = [i_w, i_t, j_w, j_t]

        # Assemble into global K
        for a in range(4):
            for b in range(4):
                K[idx[a], idx[b]] += ke[a, b]

    # ---- Build global load vector f (ndof) ----
    f = np.zeros(ndof, dtype=float)

    # Helper: collapse a point load/moment at x_p to the two nearest nodes
    def add_nodal_load(x_p: float, Fw: float, Mtheta: float) -> None:
        idx = np.searchsorted(x_nodes, x_p) - 1
        if idx < 0:
            # Left of first node -> put all on first node
            node = 0
            f[dof_per_node * node]     += Fw
            f[dof_per_node * node + 1] += Mtheta
        elif idx >= n_nodes - 1:
            # Right of last node -> put all on last node
            node = n_nodes - 1
            f[dof_per_node * node]     += Fw
            f[dof_per_node * node + 1] += Mtheta
        else:
            # Between node idx and idx+1 → distribute by distance
            xL = x_nodes[idx]
            xR = x_nodes[idx + 1]
            t = (x_p - xL) / (xR - xL)
            wL = 1.0 - t
            wR = t
            # shear contribution
            f[dof_per_node * idx]         += wL * Fw
            f[dof_per_node * (idx + 1)]   += wR * Fw
            # moment contribution
            f[dof_per_node * idx + 1]     += wL * Mtheta
            f[dof_per_node * (idx + 1)+1] += wR * Mtheta

    # ---- Add shear forces and bending moments in this bending plane ----
    # Shear forces (Fy or Fz)
    for x_p, Fw in shear_pairs:
        if Fw != 0.0:
            add_nodal_load(float(x_p), float(Fw), 0.0)

    # Bending moments (My or Mz)
    for x_p, Mth in moment_pairs:
        if Mth != 0.0:
            add_nodal_load(float(x_p), 0.0, float(Mth))

    # ---- Apply boundary conditions and solve for DOFs ----
    fixed_dofs: list[int] = []
    for i, s in enumerate(supports_sorted):
        # translational DOF (w)
        if s.type[trans_bit] == "1":
            fixed_dofs.append(dof_per_node * i)
        # rotational DOF (theta)
        if s.type[rot_bit] == "1":
            fixed_dofs.append(dof_per_node * i + 1)

    all_dofs = np.arange(ndof, dtype=int)
    fixed_dofs = np.array(sorted(set(fixed_dofs)), dtype=int)
    free_dofs = np.array([d for d in all_dofs if d not in fixed_dofs], dtype=int)

    # Displacement vector
    d = np.zeros(ndof, dtype=float)

    if free_dofs.size > 0:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        f_f = f[free_dofs]
        d_f = np.linalg.solve(K_ff, f_f)
        d[free_dofs] = d_f
    else:
        d[:] = 0.0

    # ---- Compute reactions: r = K d - f ----
    q = K @ d
    r = q - f

    # ---- Map reactions back to supports ----
    for i, s in enumerate(supports_sorted):
        w_dof = dof_per_node * i
        t_dof = w_dof + 1
        s.reactions[shear_key]  = float(r[w_dof])
        s.reactions[moment_key] = float(r[t_dof])


# -------------------------------------------------
# Loaded beam wrapper
# -------------------------------------------------
@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase

    def analyze(self, points: int = 100):
        # 1) Turn all distributed forces into point forces
        self.loads.convert_dforces_to_pforces()

        # 2) Reactions in each direction
        solve_x_forces(self.beam.supports, self.loads)
        solve_transverse_reactions(self.beam, self.loads, axis="y")
        solve_transverse_reactions(self.beam, self.loads, axis="z")

        # TODO:
        all_loads = get_all_loads(self.loads, self.beam)
        # - build internal shear/moment/deflection/torsion diagrams
        #   from these reactions + LoadCase F*/M* resultants
