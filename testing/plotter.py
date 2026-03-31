import numpy as np
from beamy.setup import LoadCase, PointForce, Beam1D, Material, Support
from sectiony import Section, Geometry, Contour, Line
from beamy.analysis import plot_beam_diagram

if __name__ == "__main__":
    # User defined triangle points (y, z)
    points = [(0,0), (2,0), (1,1), (0,0)]
    segments = [Line(start=points[i], end=points[i+1]) for i in range(len(points)-1)]
    contour = Contour(segments=segments)
    geom = Geometry(contours=[contour])

    sec = Section(name="Test", geometry=geom)
    mat = Material(name="Test", E=1, G=1)
    beam = Beam1D(
        L=5.0,
        material=mat,
        section=sec,
        supports=[
            Support(x=0.0, type="111111"),
            Support(x=5.0, type="111111")
        ]
    )

    # Create a test load case
    lc = LoadCase(name="Test Case")
    lc.add_point_force(PointForce(point=np.array([2.5, 0.5, 0.2]), force=np.array([0, -1000, 500])))
    lc.add_point_force(PointForce(point=np.array([5.0, 1.0, 0.0]), force=np.array([-100, 0, 0])))

    # Extrude along X axis with loads
    plot_beam_diagram(beam, load_case=lc)
