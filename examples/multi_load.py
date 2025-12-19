"""
Multiple Load Types Example

Demonstrates a beam with point forces, distributed forces, and moments.
Shows comprehensive loading scenarios.
"""

import numpy as np
from pathlib import Path
from sectiony.library import i_section
from beamy import (
    Beam1D, Material, Support, LoadCase, 
    PointForce, DistributedForce, Moment, LoadedMember, plot_beam_diagram
)

# Create gallery directory
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.3, b=0.15, tf=0.015, tw=0.008, r=0.01)

# 2. Create Beam (6m long, multiple supports)
L = 6.0
beam = Beam1D(
    L=L,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),   # Pinned
        Support(x=3.0, type="011000"),   # Roller
        Support(x=L, type="111000"),     # Pinned
    ]
)

# 3. Apply Multiple Load Types
loads = LoadCase(name="Multiple Loads")

# Point force at 1.5m
loads.add_point_force(PointForce(
    point=np.array([1.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -8_000.0])  # 8 kN downward
))

# Distributed force from 4m to 6m
loads.add_distributed_force(DistributedForce(
    start_position=np.array([4.0, 0.0, 0.0]),
    end_position=np.array([6.0, 0.0, 0.0]),
    start_force=np.array([0.0, 0.0, -2_000.0]),
    end_force=np.array([0.0, 0.0, -4_000.0])
))

# Moment at 2m
loads.add_moment(Moment(
    x=2.0,
    moment=np.array([0.0, 0.0, 5_000.0])  # 5 kN⋅m
))

# 4. Solve
lb = LoadedMember(beam, loads)

# 5. Get Results
print("Multiple Load Types - Results:")
print(f"Max Deflection: {lb.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb.bending('z').action.abs_max:.2f} N⋅m")
print(f"Max Shear Force: {lb.shear('z').action.abs_max:.2f} N")
print(f"Max Von Mises Stress: {lb.von_mises().max/1e6:.2f} MPa")

# 6. Visualize
plot_beam_diagram(
    lb,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "multi_load.svg")
)

print(f"\nPlot saved to: {gallery_dir / 'multi_load.svg'}")

