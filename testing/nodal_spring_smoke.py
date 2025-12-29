"""Smoke test: nodal spring stiffness affects displacements and design-grade reactions.

This is a lightweight script (no pytest) intended to be run directly:
    python testing/nodal_spring_smoke.py
"""

import numpy as np

from beamy import LoadCase, Material, NodalSpring
from beamy import Frame, Member, Node
from sectiony.library import i as i_section


def main() -> None:
    mat = Material(name="Steel", E=200e9, G=80e9)
    sect = i_section(d=0.5, b=0.2, tf=0.02, tw=0.01, r=0.0)

    L = 4.0
    members = [Member("M1", 
                     start=np.array([0.0, 0.0, 0.0]),
                     end=np.array([L, 0.0, 0.0]),
                     section=sect,
                     material=mat,
                     orientation=np.array([0.0, 1.0, 0.0]))]
    frame = Frame.from_members(members)
    
    # Apply fixed support at start node (auto-generated as N0)
    frame.nodes["N0"].support = "111111"

    P = 10_000.0  # N downward
    k_s = 2.0e6   # N/m vertical spring at B (Uz)

    loads = LoadCase("tip load + spring")
    # End node is auto-generated as N1
    loads.add_nodal_force("N1", np.array([0.0, 0.0, -P]), coords="global")

    K6 = np.zeros((6, 6))
    K6[2, 2] = k_s
    loads.nodal_springs.append(NodalSpring(node_id="N1", K=K6, coords="global"))

    frame.analyze(loads)

    dB = frame.analysis_result.nodal_displacements["N1"]
    w = float(dB[2])  # Uz

    # For a cantilever with end load in Uz, effective beam spring stiffness is k_beam = 3 E Iy / L^3
    k_beam = 3.0 * float(mat.E) * float(sect.Iy) / (L**3)
    w_expected = -P / (k_beam + k_s)

    # Reactions: beamy reports reactions as the force the structure applies to the support.
    # With a downward applied load (negative Uz), the reported Uz reaction is also negative.
    RAz = float(frame.analysis_result.reactions_physical["N0"][2])
    RAz_expected = -P * (k_beam / (k_beam + k_s))

    assert np.isfinite(w)
    assert abs(w - w_expected) <= max(1e-9, 1e-3 * abs(w_expected))
    assert abs(RAz - RAz_expected) <= max(1e-6, 1e-3 * abs(RAz_expected))

    # Sanity: equilibrium at global level (forces on the structure)
    spring_force = -k_s * w  # spring-on-structure, upward positive
    applied = -P             # applied nodal force at B in Uz
    support_on_structure = -RAz
    assert abs(applied + support_on_structure + spring_force) <= 1e-3 * P

    print("OK")
    print("  w_Bz      =", w)
    print("  w_expected=", w_expected)
    print("  R_Az      =", RAz)
    print("  R_Az_exp  =", RAz_expected)
    print("  spring_Fz =", spring_force)


if __name__ == "__main__":
    main()
