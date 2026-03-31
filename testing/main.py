# main.py
import numpy as np
from beamy.setup import Beam1D, Material, Support, LoadCase, PointForce
from sectiony import Section, Geometry, Contour, Line
from beamy.analysis import LoadedBeam

def main():
    # 1. Define Properties
    mat = Material(name="Steel", E=200e9, G=80e9)

    # Define geometry for plotting (Rectangular 0.2 x 0.1)
    # Centered at (0,0) so +/- 0.05 and +/- 0.1
    # y is vertical (height 0.2), z is horizontal (width 0.1)
    points = [
        (0.1, 0.05),   # Top Right
        (-0.1, 0.05),  # Bottom Right
        (-0.1, -0.05), # Bottom Left
        (0.1, -0.05),  # Top Left
        (0.1, 0.05)    # Close loop
    ]
    segments = [Line(start=points[i], end=points[i+1]) for i in range(len(points)-1)]
    contour = Contour(segments=segments)
    geom = Geometry(contours=[contour])

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
    P = -10_000.0

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
