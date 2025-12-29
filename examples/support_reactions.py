"""
Support Reactions Example (frame-first API)

Demonstrates how to access reactions for a single member with multiple supports
using member point supports (auto-splitting inside the solver).
"""

from __future__ import annotations

import numpy as np
from sectiony.library import i as i_section

from beamy import LoadCase, LoadedMember, Material, MemberPointForce, MemberPointMoment, MemberPointSupport


def _node_id_at_x(chain: list[tuple[float, str]], x_target: float, tol: float = 1e-6) -> str:
    for x, nid in chain:
        if abs(float(x) - float(x_target)) <= tol:
            return nid
    raise ValueError(f"No solver node found at x={x_target} (tol={tol})")


def main() -> None:
    steel = Material(name="Steel", E=200e9, G=80e9)
    section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

    L = 10.0

    # Loads: chosen to generate reactions in all directions
    lc = LoadCase("Multi-Axial Loads")
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=2.5,
            force=np.array([0.0, 0.0, -10_000.0]),  # global -Z
            coords="global",
            position_type="absolute",
        )
    )
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=7.5,
            force=np.array([0.0, 5_000.0, 0.0]),  # global +Y
            coords="global",
            position_type="absolute",
        )
    )
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=2.5,
            force=np.array([2_000.0, 0.0, 0.0]),  # global +X (axial)
            coords="global",
            position_type="absolute",
        )
    )
    lc.member_point_moments.append(
        MemberPointMoment(
            member_id="M1",
            position=7.5,
            moment=np.array([1_000.0, 0.0, 0.0]),  # global +X torsion
            coords="global",
            position_type="absolute",
        )
    )

    # Supports at x=0, x=5, x=10
    lb = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111111",  # fixed
        support_end="111000",  # pinned-like
        point_supports=[MemberPointSupport(position=5.0, support="011000")],
        load_case=lc,
    )

    frame = lb.frame
    solve_state = frame._solve_state
    if solve_state is None:
        raise RuntimeError("Expected LoadedMember to analyze on construction.")

    chain = solve_state.member_nodes_along["M1"]
    nid_0 = _node_id_at_x(chain, 0.0)
    nid_5 = _node_id_at_x(chain, 5.0)
    nid_10 = _node_id_at_x(chain, 10.0)

    result = lb.analysis_result

    def reaction_vec(nid: str) -> np.ndarray:
        if nid in result.reactions:
            return result.reactions[nid]
        return np.zeros(6, dtype=float)

    supports = [(0.0, nid_0), (5.0, nid_5), (10.0, nid_10)]

    print("\n--- Support Reactions ---\n")
    print(
        f"{'Position (m)':<15} {'Fx (N)':<12} {'Fy (N)':<12} {'Fz (N)':<12} "
        f"{'Mx (N·m)':<12} {'My (N·m)':<12} {'Mz (N·m)':<12}"
    )
    print("-" * 96)

    for x, nid in supports:
        r = reaction_vec(nid)
        print(f"{x:<15.2f} {r[0]:<12.1f} {r[1]:<12.1f} {r[2]:<12.1f} {r[3]:<12.1f} {r[4]:<12.1f} {r[5]:<12.1f}")

    print("\nNote: Reaction signs follow the global coordinate system directions.")


if __name__ == "__main__":
    main()

