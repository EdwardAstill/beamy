"""
Simple Simply Supported Beam Example

A basic example showing a simply supported beam with a point load at mid-span.
Demonstrates basic beam setup, loading, and visualization.
"""

import numpy as np
from pathlib import Path
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedMember, plot_beam_diagram

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Create gallery directory for this category
gallery_dir = PROJECT_ROOT / "gallery" / "basics"
gallery_dir.mkdir(parents=True, exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# 2. Create Beam (5m long, simply supported)
L = 5.0
beam = Beam1D(
    L=L,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111100"),  # Pinned
        Support(x=L, type="011000")     # Roller
    ]
)

# 3. Apply Loads
loads = LoadCase(name="Point Load at Mid-Span")
loads.add_point_force(PointForce(
    point=np.array([L/2, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])  # 10 kN downward
))

# 4. Solve
lb = LoadedMember(beam, loads)

# 5. Get Results
print(f"Max Deflection: {lb.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb.bending('z').action.abs_max:.2f} N⋅m")
print(f"Max Shear Force: {lb.shear('z').action.abs_max:.2f} N")

# 6. Visualize
plot_beam_diagram(
    lb,
    plot_stress=False,
    plot_section=True,
    save_path=str(gallery_dir / "simple_beam.svg")
)

print(f"\nPlot saved to: {gallery_dir / 'simple_beam.svg'}")