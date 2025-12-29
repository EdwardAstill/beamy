"""
AISC Chapter F Check Example

Demonstrates AISC 9th Edition (ASD) bending and shear checks with automatic
unit conversion. The beam is defined in SI units (meters, Pascals, Newtons),
and the check converts to AISC units (inches, ksi, kips) internally.

This example shows:
- Creating a beam with material yield strength (Fy)
- Loading the beam with vertical loads
- Running AISC Chapter F checks with unit conversion
- Interpreting the results
"""

from __future__ import annotations

import numpy as np
from sectiony.library import i as i_section

from beamy import LoadCase, LoadedMember, Material, MemberDistributedForce, MemberPointForce
from beamy.checks import aisc_9


def main() -> None:
    steel = Material(
        name="ASTM A992",
        E=200e9,
        G=80e9,
        Fy=250e6
    )

    section = i_section(
        d=0.204,
        b=0.133,
        tf=0.013,
        tw=0.008,
        r=0.010
    )

    L = 6.0

    loads = LoadCase(name="Service Loads")

    # Distributed load (global -Z) on the first half of the member
    loads.member_distributed_forces.append(
        MemberDistributedForce(
            member_id="M1",
            start_position=0.0,
            end_position=L / 2.0,
            start_force=np.array([0.0, 0.0, -5_000.0]),
            end_force=np.array([0.0, 0.0, -5_000.0]),
            coords="global",
        )
    )

    # Point load at 3L/4 (global -Z)
    loads.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=3.0 * L / 4.0,
            force=np.array([0.0, 0.0, -15_000.0]),
            coords="global",
            position_type="absolute",
        )
    )

    loaded_beam = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111100",
        support_end="011000",
        load_case=loads,
    )

    profile = loaded_beam.member_demand().actions(points=801)
    results = aisc_9.aisc_9_check(profile, length_unit="m", force_unit="N")

    print(f"AISC 9th Edition Check Results for {section.name}:")
    print(f"  Overall Pass: {results.pass_}")
    print(f"  Overall Max Utilisation: {results.utilisation:.3f}")

    # You can still get the dictionary representation if needed:
    # results_dict = results.info()


if __name__ == "__main__":
    main()
