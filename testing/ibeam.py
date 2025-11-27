import numpy as np
from sectiony.library import i_section
from beamy import Beam1D, Material, Support, LoadCase, PointForce, Moment, LoadedBeam
from beamy.analysis.plotter import plot_beam_diagram

# Create an I-beam section
# d=200mm depth, b=100mm width, tf=10mm flange, tw=6mm web, r=8mm fillet
section = i_section(d=200, b=100, tf=10, tw=6, r=8)

# Define material (steel)
steel = Material(name="Steel", E=210e9, G=80e9)

# Define supports (pinned at start, roller at end)
supports = [
    Support(x=0, type="111100"),      # Pinned (translations fixed, rotations free about y,z)
    Support(x=3000, type="011100"),   # Roller (y,z translations fixed)
]

# Create the beam (3m long)
beam = Beam1D(L=3000, material=steel, section=section, supports=supports)

# Create load case with point forces
loads = LoadCase(name="Test Load")
loads.add_point_force(PointForce(
    point=np.array([1500, 0, 0]),   # At midspan
    force=np.array([0, -1000, 1000])  # 10kN downward
))
loads.add_point_force(PointForce(
    point=np.array([1000, 0, 0]),   # At 1m
    force=np.array([0, -5000, 0])   # 5kN downward
))
loads.add_moment(Moment(
    x=0,
    moment=np.array([0, 10000, 1000])
))

# Create loaded beam (solves for reactions)
lb = LoadedBeam(beam, loads)

# Plot the beam with forces and stress coloring
plot_beam_diagram(lb, plot_stress=True, plot_section=True)
