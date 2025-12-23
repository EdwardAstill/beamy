import numpy as np
from pathlib import Path
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, Moment, LoadedMember, DistributedForce, plot_beam_diagram

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Create an I-beam section
# d=200mm depth, b=100mm width, tf=10mm flange, tw=6mm web, r=8mm fillet
section = i_section(d=200, b=100, tf=10, tw=6, r=8)

# Define material (steel)
steel = Material(name="Steel", E=210e9, G=80e9)

# Define supports (pinned at start, roller at end)
supports = [
    Support(x=0, type="111100"),      # Pinned (translations fixed, rotations free about y,z)
    Support(x=1300, type="111111"),  # Pinned (translations fixed, rotations free about y,z)
    Support(x=3000, type="011100"),   # Roller (y,z translations fixed)
]

# Create the beam (3m long)
beam = Beam1D(L=3000, material=steel, section=section, supports=supports)

# Create load case with point forces
loads = LoadCase(name="Test Load")

loads.add_point_force(PointForce(
    point=np.array([1000, 0, 0]),   # At 1m
    force=np.array([0, -5000, 0])   # 5kN downward
))
# loads.add_moment(Moment(
#     x=0,
#     moment=np.array([0, 10000, 1000])
# ))
loads.add_distributed_force(DistributedForce(
    start_position=np.array([2500, 0, 0]),
    end_position=np.array([1300, 0, 0]),
    start_force=np.array([0, -1, 0]),
    end_force=np.array([0, -5, 0])
))

# Create loaded beam (solves for reactions)
lb = LoadedMember(beam, loads)

# Find the location of maximum Von Mises stress
vm_results = lb.von_mises(points=200)  # Use more points for better resolution
max_vm_idx = np.argmax(vm_results._values)
max_vm_x = vm_results._x[max_vm_idx]
max_vm_val = vm_results._values[max_vm_idx]

print(f"Max Von Mises Stress: {max_vm_val/1e6:.2f} MPa at x = {max_vm_x:.2f} mm")

# Plot the beam with forces and stress coloring
from beamy import StressPlotter

gallery_dir = PROJECT_ROOT / "gallery" / "section"
gallery_dir.mkdir(parents=True, exist_ok=True)

# Save beam diagram
plot_beam_diagram(
    lb,
    plot_stress=True,
    plot_section=True,
    save_path=str(gallery_dir / "ibeam_3d_diagram.svg")
)

# Plot section stress at the critical location
sp = StressPlotter(lb)
sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="von_mises",
    title=f"Max Von Mises Stress Section (x={max_vm_x:.0f}mm)",
    cmap="plasma",
    show=False
)
import matplotlib.pyplot as plt
plt.savefig(gallery_dir / "ibeam_section_stress.svg", bbox_inches='tight', dpi=300)
plt.close()