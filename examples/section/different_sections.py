"""
Different Section Types Example

Demonstrates beams with different cross-section shapes:
- I-beam
- Rectangular section
- Circular section (if available)
Shows how different sections affect beam behavior.
"""

import numpy as np
from pathlib import Path
from sectiony.library import i_section, rhs
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedMember, plot_beam_diagram

PROJECT_ROOT = Path(__file__).resolve().parents[2]

gallery_dir = PROJECT_ROOT / "gallery" / "section"
gallery_dir.mkdir(parents=True, exist_ok=True)

steel = Material(name="Steel", E=200e9, G=80e9)
L = 4.0

# Common loading
loads = LoadCase(name="Point Load")
loads.add_point_force(PointForce(
    point=np.array([L/2, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])  # 10 kN downward
))

# Common supports
supports = [
    Support(x=0.0, type="111100"),
    Support(x=L, type="011010")
]

# 1. I-Beam
print("=== I-Beam ===")
section_ibeam = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.008)
beam_ibeam = Beam1D(L=L, material=steel, section=section_ibeam, supports=supports)
lb_ibeam = LoadedMember(beam_ibeam, loads)

print(f"Max Deflection: {lb_ibeam.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb_ibeam.bending('z').action.abs_max:.2f} N⋅m")

plot_beam_diagram(
    lb_ibeam,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "section_ibeam.svg")
)

# 2. Rectangular Hollow Section
print("\n=== Rectangular Hollow Section ===")
section_rect = rhs(b=0.1, h=0.2, t=0.005, r=0.0)  # 100mm x 200mm RHS with 5mm wall
beam_rect = Beam1D(L=L, material=steel, section=section_rect, supports=supports)
lb_rect = LoadedMember(beam_rect, loads)

print(f"Max Deflection: {lb_rect.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb_rect.bending('z').action.abs_max:.2f} N⋅m")

plot_beam_diagram(
    lb_rect,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "section_rhs.svg")
)

print(f"\nPlots saved to gallery directory:")
print(f"  - {gallery_dir / 'section_ibeam.svg'}")
print(f"  - {gallery_dir / 'section_rhs.svg'}")