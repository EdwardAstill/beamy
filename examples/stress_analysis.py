"""
Stress Analysis Example

Demonstrates comprehensive stress analysis including:
- Finding critical stress locations
- Plotting section stresses at multiple locations
- Comparing different stress types
"""

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from sectiony.library import i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam, StressPlotter
from beamy.analysis.beam_plotter import plot_beam_diagram

# Create gallery directory
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)

# 1. Define Properties
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.25, b=0.12, tf=0.012, tw=0.007, r=0.01)

# 2. Create Beam
L = 4.0
beam = Beam1D(
    L=L,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),
        Support(x=L, type="011000")
    ]
)

# 3. Apply Loads
loads = LoadCase(name="Stress Analysis")
loads.add_point_force(PointForce(
    point=np.array([1.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -12_000.0])  # 12 kN downward
))
loads.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -8_000.0])   # 8 kN downward
))

# 4. Solve
lb = LoadedBeam(beam, loads)

# 5. Find Critical Stress Location
vm_results = lb.von_mises(points=200)
max_vm_idx = np.argmax(vm_results._values)
max_vm_x = vm_results._x[max_vm_idx]
max_vm_val = vm_results._values[max_vm_idx]

print("Stress Analysis Results:")
print(f"Max Von Mises Stress: {max_vm_val/1e6:.2f} MPa")
print(f"Location of Max Stress: x = {max_vm_x:.3f} m")

# 6. Plot 3D Beam Diagram with Stress
plot_beam_diagram(
    lb,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "stress_analysis_3d.svg")
)

# 7. Plot Section Stresses
sp = StressPlotter(lb)

# Plot at maximum stress location
sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="von_mises",
    title=f"Von Mises Stress at Critical Location (x={max_vm_x:.3f}m)",
    cmap="plasma",
    show=False
)
plt.savefig(gallery_dir / "stress_analysis_von_mises.svg", bbox_inches='tight', dpi=300)
plt.close()

# Plot bending stress at mid-span
sp.plot_stress_at(
    x_pos=L/2,
    stress_type="sigma_bending",
    title=f"Bending Stress at Mid-Span (x={L/2:.1f}m)",
    cmap="RdBu_r",
    show=False
)
plt.savefig(gallery_dir / "stress_analysis_bending.svg", bbox_inches='tight', dpi=300)
plt.close()

# Plot multiple stress types at same location
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="sigma_bending",
    ax=axes[0],
    show=False,
    cmap="RdBu_r"
)
axes[0].set_title("Bending Stress")

sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="tau_shear",
    ax=axes[1],
    show=False,
    cmap="viridis"
)
axes[1].set_title("Shear Stress")

sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="von_mises",
    ax=axes[2],
    show=False,
    cmap="plasma"
)
axes[2].set_title("Von Mises Stress")

plt.suptitle(f"Stress Components at x={max_vm_x:.3f}m", fontsize=14)
plt.tight_layout()
plt.savefig(gallery_dir / "stress_analysis_components.svg", bbox_inches='tight', dpi=300)
plt.close()

print(f"\nPlots saved to gallery directory:")
print(f"  - stress_analysis_3d.svg")
print(f"  - stress_analysis_von_mises.svg")
print(f"  - stress_analysis_bending.svg")
print(f"  - stress_analysis_components.svg")

