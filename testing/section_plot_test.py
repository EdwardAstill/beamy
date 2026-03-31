import matplotlib.pyplot as plt
from sectiony import Section, Geometry, Contour, Line

def _make_contour(points: list[tuple[float, float]], close: bool = True) -> Contour:
    """Helper to create a Contour from a list of points."""
    if close and points[0] != points[-1]:
        points = list(points) + [points[0]]
    segments = [Line(start=points[i], end=points[i+1]) for i in range(len(points)-1)]
    return Contour(segments=segments)

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

    contour_outer = _make_contour(points_outer)
    contour_inner = _make_contour(points_inner)

    geom = Geometry(contours=[contour_outer, contour_inner])
    sec = Section(name="Hollow Rect", geometry=geom)

    # Plot showing the plot
    fig, ax = plt.subplots()
    sec.plot(ax=ax, show=False)

    print("  [PASS] Section plot executed")

if __name__ == "__main__":
    test_section_plot()
