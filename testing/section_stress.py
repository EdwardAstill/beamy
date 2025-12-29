import sys
from pathlib import Path

# Add the src directory to the path so we can import beamy as a package
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import numpy as np
import matplotlib.pyplot as plt
from sectiony.library import i_section
from sectiony.stress import Stress
from beamy import LoadCase, LoadedMember, Material, MemberPointForce, MemberPointMoment

def test_section_stress():
    # 1. Setup Beam
    # I-beam: d=200, b=100, tf=10, tw=6
    section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)
    steel = Material(name="Steel", E=210e9, G=80e9)
    
    L = 3.0
    
    # 2. Loads: Combined bending, axial, and torsion
    loads = LoadCase(name="Combined")
    loads.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=1.5,
            force=np.array([10000.0, -5000.0, 2000.0]),
            coords="global",
            position_type="absolute",
        )
    )
    loads.member_point_moments.append(
        MemberPointMoment(
            member_id="M1",
            position=1.5,
            moment=np.array([1000.0, 5000.0, -2000.0]),
            coords="global",
            position_type="absolute",
        )
    )
    
    lb = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111100",
        support_end="011100",
        load_case=loads,
    )
    
    # 3. Analyze at specific location x_loc
    x_loc = 1.4 # Just before the load
    
    print(f"Analyzing section stress at x = {x_loc}")
    
    # Get internal actions at x_loc
    # Note: .at(x) interpolates. For step changes (point loads), be careful with exact location.
    N  = lb.axial().action.at(x_loc)      # Axial force
    Vy = lb.shear('y').action.at(x_loc)   # Shear force in Y
    Vz = lb.shear('z').action.at(x_loc)   # Shear force in Z
    Mx = lb.torsion().action.at(x_loc)    # Torsional moment
    My = lb.bending('y').action.at(x_loc) # Bending moment about Y
    Mz = lb.bending('z').action.at(x_loc) # Bending moment about Z
    
    print(f"Internal Actions:")
    print(f"  N  (Axial) = {N:.2f}")
    print(f"  Vy (Shear Y) = {Vy:.2f}")
    print(f"  Vz (Shear Z) = {Vz:.2f}")
    print(f"  Mx (Torque) = {Mx:.2f}")
    print(f"  My (Bend Y) = {My:.2f}")
    print(f"  Mz (Bend Z) = {Mz:.2f}")
    
    # 4. Create Stress object using sectiony
    stress = Stress(
        section=section,
        N=N,
        Vy=Vy,
        Vz=Vz,
        Mx=Mx,
        My=My,
        Mz=Mz
    )
    
    # 5. Plot using sectiony's built-in plotting
    # Note: sectiony plots with Y horizontal, Z vertical by default
    # We want Z horizontal, Y vertical, so we'll need to handle the axes
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Use sectiony's plot method
    stress.plot(
        stress_type="total",
        ax=ax,
        show=False,
        cmap="plasma"
    )
    
    # Adjust axes labels and orientation
    # sectiony plots (y, z) -> (x, y) in matplotlib
    # We want z on horizontal, y on vertical
    # So we need to swap the axes
    ax.set_xlabel("Z")
    ax.set_ylabel("Y")
    ax.set_title(f"Section Stress at x={x_loc}\n(Axial + Bending + Torsion)")
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    test_section_stress()
