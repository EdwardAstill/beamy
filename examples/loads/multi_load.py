"""
Multiple Load Types Example

Demonstrates a beam with point forces, distributed forces, and moments.
Shows comprehensive loading scenarios.
"""

import numpy as np
from pathlib import Path
from sectiony.library import i_section
from beamy import (
    Material,
    LoadCase,
    MemberPointForce,
    MemberPointMoment,
    MemberDistributedForce,
    MemberPointSupport,
    LoadedMember,
    plot_beam_diagram,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Create gallery directory for this category
gallery_dir = PROJECT_ROOT / "gallery" / "loads"
gallery_dir.mkdir(parents=True, exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.3, b=0.15, tf=0.015, tw=0.008, r=0.01)

# 2. Create member (6m long) + supports (supports along the member belong to geometry)
L = 6.0

# 3. Apply Multiple Load Types
load_case = LoadCase(name="Multiple Loads")

# Point force at 1.5m
load_case.member_point_forces.append(
    MemberPointForce(
        member_id="M1",
        position=1.5,
        force=np.array([0.0, 0.0, -8_000.0]),  # 8 kN downward
        coords="global",
        position_type="absolute",
    )
)

# Distributed force from 4m to 6m
load_case.member_distributed_forces.append(
    MemberDistributedForce(
        member_id="M1",
        start_position=4.0,
        end_position=6.0,
        start_force=np.array([0.0, 0.0, -2_000.0]),
        end_force=np.array([0.0, 0.0, -4_000.0]),
        coords="global",
    )
)

# Moment at 2m
load_case.member_point_moments.append(
    MemberPointMoment(
        member_id="M1",
        position=2.0,
        moment=np.array([0.0, 0.0, 5_000.0]),  # 5 kN⋅m about global z
        coords="global",
        position_type="absolute",
    )
)

# 4. Solve
lb = LoadedMember(
    id="M1",
    start=np.array([0.0, 0.0, 0.0]),
    end=np.array([L, 0.0, 0.0]),
    section=section,
    material=steel,
    orientation=np.array([0.0, 1.0, 0.0]),  # local z aligns with global z
    support_start="111000",
    support_end="111000",
    point_supports=[MemberPointSupport(position=3.0, support="011000")],
    load_case=load_case,
)

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