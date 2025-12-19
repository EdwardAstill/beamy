import numpy as np
from sectiony.library import i as i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, Moment, LoadedMember, DistributedForce
from beamy.analysis import plot_analysis_results
from pathlib import Path

# 1. Setup Beam
# ---------------------------------------------
# Create an I-beam section
section = i_section(d=200, b=100, tf=10, tw=6, r=8)

# Define material (steel)
steel = Material(name="Steel", E=210e9, G=80e9)

# Define supports
# Pinned at start (x=0), Roller at end (x=3000)
# Note: Internal nodes are automatically inserted at load positions for accurate deflection
supports = [
    Support(x=0, type="111100"),      # Pinned (Fixed: Ux, Uy, Uz, Rx)
    Support(x=3000, type="011100"),   # Roller (Fixed: Uy, Uz, Rx)
]

# Create the beam (3m long)
beam = Beam1D(L=3000, material=steel, section=section, supports=supports)

# 2. Setup Loads
# ---------------------------------------------
loads = LoadCase(name="Complex Loading")

# Vertical Point Force (Fy - Bending about Z)
loads.add_point_force(PointForce(
    point=np.array([1000, 0, 0]),   # At x=1m
    force=np.array([0, -5000, 0])   # 5kN downward
))

# Horizontal Point Force (Fz - Bending about Y)
loads.add_point_force(PointForce(
    point=np.array([2000, 0, 0]),   # At x=2m
    force=np.array([0, 0, 2000])    # 2kN sideways
))

# Distributed Force (Vertical)
loads.add_distributed_force(DistributedForce(
    start_position=np.array([1500, 0, 0]),
    end_position=np.array([2500, 0, 0]),
    start_force=np.array([0, -2000, 0]), # -2 kN/m
    end_force=np.array([0, -2000, 0])    # -2 kN/m
))

# Axial Force (Compression)
loads.add_point_force(PointForce(
    point=np.array([3000, 0, 0]),   # At end
    force=np.array([-10000, 0, 0])  # 10kN axial compression
))

# Torsion
loads.add_moment(Moment(
    x=1500,
    moment=np.array([500, 0, 0])    # 500 Nm Torsion
))

# 3. Analyze and Plot
# ---------------------------------------------
lb = LoadedMember(beam, loads)

# Create gallery directory if it doesn't exist
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)

print("Generating analysis plots...")
plot_analysis_results(
    lb, 
    save_path=str(gallery_dir / "line_plots_example.svg"),
    show=False,  # Set to True to display interactively
    units={'length': 'mm', 'force': 'N', 'moment': 'N.mm', 'deflection': 'mm'}
)
print(f"Plots saved to {gallery_dir / 'line_plots_example.svg'}")

