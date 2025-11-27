# main.py
import sys
from pathlib import Path

# Add the src directory to the path so we can import beamy as a package
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import numpy as np
from beamy.setup import Beam1D, Material, Support, LoadCase, PointForce
from beamy.section import Section, Geometry, Shape
from beamy.analysis import LoadedBeam

def main():
    # 1. Define Properties
    mat = Material(name="Steel", E=200e9, G=80e9)
    
    # Define geometry for plotting (Rectangular 0.2 x 0.1)
    # Centered at (0,0) so +/- 0.05 and +/- 0.1
    # y is vertical (height 0.2), z is horizontal (width 0.1)
    # y range: -0.1 to 0.1
    # z range: -0.05 to 0.05
    points = [
        (0.1, 0.05),   # Top Right
        (-0.1, 0.05),  # Bottom Right
        (-0.1, -0.05), # Bottom Left
        (0.1, -0.05),  # Top Left
        (0.1, 0.05)    # Close loop
    ]
    shape = Shape(points=points)
    geom = Geometry(shapes=[shape])
    
    sec = Section(name="Test", geometry=geom) # Properties calculated from geometry
    
    L = 1.0  # 1 m beam

    # 2. Create Beam
    beam = Beam1D(
        L=L,
        material=mat,
        section=sec,
        supports=[
            Support(x=0.0, type="111100"), # Pinned + Fixed Rotation about X
            Support(x=L, type="011000")    # Roller
        ]
    )

    # 3. Apply Loads
    # 10 kN downward at midspan
    P = -10_000.0  # negative in local z (since we plot z as horizontal transverse usually, but let's see. 
    # beam.py: y is usually vertical bending axis (Iz), z is horizontal (Iy). 
    # Force in z is horizontal load. Force in y is vertical load.
    
    load = PointForce(
        point=np.array([0.5, 0.0, 0.0]), # x, y, z
        force=np.array([0.0, 0.0, P]), # Load in Z direction (Transverse horizontal)
    )

    lc = LoadCase(name="Case 1")
    lc.add_point_force(load)

    # 4. Solve
    lb = LoadedBeam(beam, lc)

    # 5. Output Results
    print("All Loads (Applied + Reactions):")
    for x, t, v in lb.all_loads:
        print(f"  x={x:.2f}, type={t}, value={v:.2f}")
        
    print("\nResults along beam:")
    shear_res = lb.shear("z")
    bend_res = lb.bending("z")
    defl_res = lb.deflection("z")
    
    for x in np.linspace(0.0, L, 11):
        vz = shear_res.action.at(x)
        mz = bend_res.action.at(x)
        w = defl_res.at(x)
        print(f"x={x:.2f}: Vz={vz:.2f}, Mz={mz:.2f}, w={w:.6e}")

    # 6. Plot
    print("\nPlotting beam diagram...")
    lb.plot()

if __name__ == "__main__":
    main()

