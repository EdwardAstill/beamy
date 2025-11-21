import sys
from pathlib import Path
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import numpy as np
from beamy.setup import LoadCase, PointForce, Beam1D, Material, Support
from beamy.section import Section, Geometry, Shape
from beamy.analysis import plot_beam_diagram

if __name__ == "__main__":
    # User defined triangle points (y, z)
    points = ((0,0), (2,0), (1,1), (0,0))
    
    # Create beam objects
    shape = Shape(points=points)
    geom = Geometry(shapes=[shape])
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
    # Add a point force at x=2.5 (mid-span), y=0.5, z=0.2
    # Force vector pointing down and sideways
    lc.add_point_force(PointForce(point=np.array([2.5, 0.5, 0.2]), force=np.array([0, -1000, 500])))
    # Add another force at end
    lc.add_point_force(PointForce(point=np.array([5.0, 1.0, 0.0]), force=np.array([-100, 0, 0])))
    
    # Extrude along X axis with loads
    plot_beam_diagram(beam, load_case=lc)
