"""
Support Reactions Example

Demonstrates how to access reaction forces and moments at supports
for a beam with multiple supports and complex loading.
"""

import numpy as np
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, Moment, LoadedBeam

# 1. Define Material and Section
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# 2. Define Beam with Multiple Supports (Continuous Beam)
# Length = 10m
# Supports at x=0, x=5, x=10
beam = Beam1D(
    L=10.0,
    material=steel,
    section=section,
    supports=[
        # Left: Fully Fixed (restrains everything)
        Support(x=0.0, type="111111"),
        
        # Middle: Roller in Z (restrains Y and Z translation, free rotation)
        # Note: Beamy support string is [Ux, Uy, Uz, Rx, Ry, Rz]
        # "011000" -> Free X, Fixed Y, Fixed Z, Free rotations
        Support(x=5.0, type="011000"),
        
        # Right: Pinned (Fixed X, Y, Z translation, free rotation)
        Support(x=10.0, type="111000")
    ]
)

# 3. Apply Complex Loads
# We apply loads that will generate reactions in all directions
lc = LoadCase("Multi-Axial Loads")

# Vertical Load (Generates Fz, My reactions)
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0]) # -10 kN in Z
))

# Horizontal Transverse Load (Generates Fy, Mz reactions)
lc.add_point_force(PointForce(
    point=np.array([7.5, 0.0, 0.0]),
    force=np.array([0.0, 5_000.0, 0.0])   # +5 kN in Y
))

# Axial Load (Generates Fx reactions)
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([2_000.0, 0.0, 0.0])   # +2 kN in X
))

# Torsion (Generates Mx reactions)
lc.add_moment(Moment(
    x=7.5,
    moment=np.array([1_000.0, 0.0, 0.0])  # 1 kNm Torque
))

# 4. Analyze
lb = LoadedBeam(beam, lc)

# 5. Access and Print Reactions
print("\n--- Support Reactions ---\n")
print(f"{'Position (m)':<15} {'Fx (N)':<12} {'Fy (N)':<12} {'Fz (N)':<12} {'Mx (Nm)':<12} {'My (Nm)':<12} {'Mz (Nm)':<12}")
print("-" * 90)

# Sort supports by position for nice output
sorted_supports = sorted(beam.supports, key=lambda s: s.x)

for s in sorted_supports:
    # Reactions are stored in the dictionary s.reactions
    rx = s.reactions
    
    print(f"{s.x:<15.2f} {rx.get('Fx', 0):<12.1f} {rx.get('Fy', 0):<12.1f} {rx.get('Fz', 0):<12.1f} {rx.get('Mx', 0):<12.1f} {rx.get('My', 0):<12.1f} {rx.get('Mz', 0):<12.1f}")

print("\nNote: Reaction signs follow the coordinate system directions.")

