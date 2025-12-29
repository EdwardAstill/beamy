"""
AISC Chapter F Check Example (Standalone)

Demonstrates AISC 9th Edition (ASD) bending and shear checks with automatic
unit conversion. The beam is defined in SI units (meters, Pascals, Newtons),
and the check converts to AISC units (inches, ksi, kips) internally.

Run from beamy root with: python -m examples.aisc_check_standalone
Or: cd beamy && python -c "import sys; sys.path.insert(0, 'src'); from examples.aisc_check_standalone import *"
"""

from __future__ import annotations

import numpy as np
from sectiony.library import i as i_section

from beamy import LoadCase, LoadedMember, Material, MemberDistributedForce, MemberPointForce
from beamy.checks import aisc_9


def main() -> None:
    steel = Material(
        name="ASTM A992",
        E=200e9,      # Young's modulus (Pa)
        G=80e9,       # Shear modulus (Pa)
        Fy=250e6      # Yield stress (Pa)
    )

    section = i_section(
        d=0.204,      # Depth: 204 mm
        b=0.133,      # Width: 133 mm
        tf=0.013,     # Flange thickness: 13 mm
        tw=0.008,     # Web thickness: 8 mm
        r=0.010       # Root radius: 10 mm
    )

    L = 6.0  # 6 meter span

    loads = LoadCase(name="Service Loads")

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


if __name__ == "__main__":
    main()

