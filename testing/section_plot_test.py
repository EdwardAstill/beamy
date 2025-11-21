import sys
from pathlib import Path
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import matplotlib.pyplot as plt
from beamy.section import Section, Geometry, Shape

def test_section_plot():
    print("Testing section plot...")
    # Define geometry (Rectangular 0.2 x 0.1) with a hole
    points_outer = [
        (0.1, 0.05),   # Top Right
        (-0.1, 0.05),  # Bottom Right
        (-0.1, -0.05), # Bottom Left
        (0.1, -0.05),  # Top Left
    ]
    
    points_inner = [
        (0.05, 0.025),
        (-0.05, 0.025),
        (-0.05, -0.025),
        (0.05, -0.025),
    ]
    
    shape_outer = Shape(points=points_outer)
    shape_inner = Shape(points=points_inner, hollow=True)
    
    geom = Geometry(shapes=[shape_outer, shape_inner])
    sec = Section(name="Hollow Rect", geometry=geom)
    
    # Plot showing the plot
    fig, ax = plt.subplots()
    sec.plot(ax=ax, show=True)
    
    print("  [PASS] Section plot executed")

if __name__ == "__main__":
    test_section_plot()

