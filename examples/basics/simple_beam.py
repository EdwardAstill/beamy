"""
Simple Simply Supported Beam Example

A basic example showing a simply supported beam with a point load at mid-span.
Demonstrates basic beam setup, loading, and visualization.
"""

import numpy as np
from pathlib import Path
from sectiony.library import i as i_section
from beamy import Material, LoadCase, MemberPointForce, LoadedMember, plot_beam_diagram

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Create gallery directory for this category
gallery_dir = PROJECT_ROOT / "gallery" / "basics"
gallery_dir.mkdir(parents=True, exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# 2. Create member (5m long, simply supported)
L = 5.0

# 3. Apply Loads
load_case = LoadCase(name="Point Load at Mid-Span")
load_case.member_point_forces.append(
    MemberPointForce(
        member_id="M1",
        position=L / 2.0,
        force=np.array([0.0, 0.0, -10_000.0]),  # 10 kN downward (global Z)
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
    support_start="111100",  # pinned-like (translations + Rx fixed)
    support_end="011000",  # roller-like (Uy/Uz fixed)
    load_case=load_case,
)

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