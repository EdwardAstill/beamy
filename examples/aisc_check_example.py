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

import numpy as np
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, DistributedForce, LoadedBeam


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

    beam = Beam1D(
        L=L,
        material=steel,
        section=section,
        supports=[
            Support(x=0.0, type="111100"),
            Support(x=L, type="011000")
        ]
    )

    loads = LoadCase(name="Service Loads")

    loads.add_distributed_force(DistributedForce(
        start_position=np.array([0.0, 0.0, 0.0]),
        end_position=np.array([L / 2, 0.0, 0.0]),
        start_force=np.array([0.0, 0.0, -5_000.0]),
        end_force=np.array([0.0, 0.0, -5_000.0])
    ))

    loads.add_point_force(PointForce(
        point=np.array([3 * L / 4, 0.0, 0.0]),
        force=np.array([0.0, 0.0, -15_000.0])
    ))

    loaded_beam = LoadedBeam(beam, loads)
    results = loaded_beam.check_aisc_chapter_f(length_unit="m", force_unit="N")

    print(results)


if __name__ == "__main__":
    main()