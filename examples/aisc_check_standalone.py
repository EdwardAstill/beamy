"""
AISC Chapter F Check Example (Standalone)

Demonstrates AISC 9th Edition (ASD) bending and shear checks with automatic
unit conversion. The beam is defined in SI units (meters, Pascals, Newtons),
and the check converts to AISC units (inches, ksi, kips) internally.

Run from beamy root with: python -m examples.aisc_check_standalone
Or: cd beamy && python -c "import sys; sys.path.insert(0, 'src'); from examples.aisc_check_standalone import *"
"""

import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, DistributedForce, LoadedBeam


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

    print(f"AISC Check Results for {loaded_beam.beam.section.name}:")
    print(f"  Overall Pass: {results.pass_}")
    print(f"  Overall Max Utilisation: {results.utilisation:.3f}")

    print("\nBending Checks:")
    for res in results.bending:
        print(f"  Axis: {res.axis}")
        print(f"  Max Utilisation: {res.utilisation:.3f}")
        for i, seg in enumerate(res.segments):
            print(f"    Segment {i+1} ({seg.x_start:.2f}m - {seg.x_end:.2f}m):")
            print(f"      Governing Fb: {seg.governing_Fb:.1f} ksi")
            print(f"      Demand Fb:    {seg.demand_Fb:.1f} ksi")
            print(f"      Utilisation:  {seg.utilisation:.3f}")

    print("\nShear Check:")
    print(f"  Max Shear (V_max): {results.shear.V_max:.1f} kips")
    if results.shear.capacity:
        print(f"  Capacity:          {results.shear.capacity:.1f} kips")
    print(f"  Utilisation:       {results.shear.utilisation:.3f}")


if __name__ == "__main__":
    main()

