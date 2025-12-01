"""
Cantilever Beam Example

A cantilever beam with a point load at the free end.
Demonstrates fixed support and end loading.
"""

import numpy as np
from pathlib import Path
from sectiony.library import rhs
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam
from beamy.analysis.beam_plotter import plot_beam_diagram

# Create gallery directory
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = rhs(b=0.1, h=0.2, t=0.005, r=0.0)  # 100mm x 200mm RHS with 5mm wall

# 2. Create Cantilever Beam (3m long, fixed at x=0)
L = 3.0
beam = Beam1D(
    L=L,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111111"),  # Fixed (all DOFs constrained)
    ]
)

# 3. Apply Load at Free End
loads = LoadCase(name="End Load")
loads.add_point_force(PointForce(
    point=np.array([L, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -5_000.0])  # 5 kN downward
))

# 4. Solve
lb = LoadedBeam(beam, loads)

# 5. Get Results
print("Cantilever Beam Results:")
print(f"Max Deflection: {lb.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb.bending('z').action.abs_max:.2f} N⋅m")
print(f"Max Shear Force: {lb.shear('z').action.abs_max:.2f} N")

# 6. Visualize
plot_beam_diagram(
    lb,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "cantilever.svg")
)

print(f"\nPlot saved to: {gallery_dir / 'cantilever.svg'}")

