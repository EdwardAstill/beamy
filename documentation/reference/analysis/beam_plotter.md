# Beam Plotter

Reference documentation for 3D visualization of beams, loads, and stress distributions.

## Table of Contents

- [plot_beam_diagram](#plot_beam_diagram)
- [plot_loads](#plot_loads)

---

## plot_beam_diagram

Plots a 3D beam diagram with section outline, loads, supports, and optional stress visualization.

```python
from beamy.analysis.beam_plotter import plot_beam_diagram
# Or import directly from beamy:
from beamy import plot_beam_diagram  # Available via __init__.py

def plot_beam_diagram(
    loaded_beam: LoadedBeam,
    plot_stress: bool = False,
    plot_section: bool = True,
    save_path: Optional[str] = None
) -> None:
    """
    Plots a 3D beam diagram with:
    - Optional 2D section outline on the YZ plane at x=0 (shear center at origin)
    - A line along the beam x-axis representing the beam length
    - Optional von Mises stress coloring on the beam axis
    - Point forces, distributed forces, and moments as 3D arrows/arcs
    - Supports as hollow circles with labels
    
    Args:
        loaded_beam: LoadedBeam object containing beam, loads, and analysis results
        plot_stress: If True, color the beam axis by von Mises stress
        plot_section: If True, draw the section outline at x=0
        save_path: If provided, save the plot to this file path instead of showing it
    """
```

### Parameters

- `loaded_beam` (LoadedBeam): LoadedBeam object containing beam, loads, and analysis results
- `plot_stress` (bool): If True, color the beam axis by von Mises stress (default: False)
- `plot_section` (bool): If True, draw the section outline at x=0 (default: True)
- `save_path` (str, optional): If provided, saves the plot to this file path instead of displaying it. Uses `bbox_inches='tight'` and `dpi=300` for high-quality output.

### Visual Features

- **Section Outline**: 2D cross-section drawn on the YZ plane at x=0, positioned with shear center at origin. The section face is filled with light grey.
- **Beam Axis**: Line along the beam length, optionally colored by von Mises stress using a plasma colormap
- **Point Forces**: Red 3D arrows with:
  - Arrow tip at the application point
  - Length proportional to force magnitude
  - Longest arrow is 1/5th of the beam length
  - 3D cone arrowhead with lighter red base
- **Distributed Forces**: Green 3D arrows showing the distributed load:
  - Multiple arrows along the distributed load length (minimum 2, scales with load length)
  - All arrow tails connected with a line
  - Fixed cone size for all arrows in the same distributed load
  - Arrows scale proportionally to force magnitude
- **Moments**: Blue 3/4 arcs with cone tips showing moment direction:
  - Arc radius proportional to moment magnitude
  - Maximum arc radius is 1/10th of the beam length
  - Cone tip indicates the direction of the moment vector
- **Supports**: Hollow black circles (markers) positioned along the beam axis with 6-digit support type labels
- **Stress Visualization**: When `plot_stress=True`, the beam axis is color-coded using a plasma colormap with a colorbar positioned close to the diagram

### Coordinate Mapping

The plot uses a specific coordinate mapping for visualization:
- **Plot X** ← Beam Z (Transverse Horizontal)
- **Plot Y** ← Beam X (Longitudinal)
- **Plot Z** ← Beam Y (Vertical)

This mapping provides a natural view where:
- The beam extends along the Y-axis (Plot Y = Beam X)
- The cross-section is viewed in the XZ plane (Plot X = Beam Z, Plot Z = Beam Y)

### Example

```python
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam
from beamy.analysis.beam_plotter import plot_beam_diagram
from sectiony.library import i_section
import numpy as np

# Create beam
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)
beam = Beam1D(
    L=5.0,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),
        Support(x=5.0, type="111000")
    ]
)

# Create loads
lc = LoadCase(name="Point Load")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Create loaded beam
loaded_beam = LoadedBeam(beam, lc)

# Plot with stress visualization
plot_beam_diagram(loaded_beam, plot_stress=True, plot_section=True)

# Save to file instead of showing
plot_beam_diagram(loaded_beam, plot_stress=True, plot_section=True, save_path="beam_diagram.svg")
```

### Notes

- The function creates a new matplotlib figure and displays it using `plt.show()` unless `save_path` is provided
- The view is set to elevation=20°, azimuth=-60° for a good default perspective
- Plot limits are automatically calculated based on section size and beam length
- The section must have geometry defined (via `section.geometry`) for section plotting to work
- All arrows and arcs use `zorder` to ensure proper layering (shafts appear in front of bases)
- The plot has a white background with no grid or panes visible

---

## plot_loads

Plots loads on an existing 3D axes object.

```python
from beamy.analysis.beam_plotter import plot_loads

def plot_loads(
    ax: matplotlib.axes.Axes3D,
    load_case: LoadCase,
    beam_length: float
) -> None:
    """
    Plots loads on a 3D axes object.
    
    Args:
        ax: 3D matplotlib axes object
        load_case: LoadCase object containing forces to plot
        beam_length: Length of the beam (used for scaling arrow lengths)
    """
```

### Parameters

- `ax` (matplotlib.axes.Axes3D): 3D axes object from `fig.add_subplot(111, projection='3d')`
- `load_case` (LoadCase): Load case containing forces to visualize
- `beam_length` (float): Length of the beam, used to scale arrow lengths

### Behavior

- Plots point forces as red 3D arrows
- Plots distributed forces as green 3D arrows with connecting lines
- Plots moments as blue 3/4 arcs with cone tips
- All loads are scaled proportionally to their magnitudes

### Arrow Scaling

- Arrow lengths are proportional to force magnitude
- The longest arrow (corresponding to the maximum force magnitude) is set to `beam_length / 5`
- All other arrows scale proportionally
- Distributed force arrows use a fixed cone size based on average force magnitude

### Example

```python
import matplotlib.pyplot as plt
from beamy.analysis.beam_plotter import plot_loads
from beamy import LoadCase, PointForce
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Create load case
lc = LoadCase(name="Forces")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Plot loads
plot_loads(ax, lc, beam_length=5.0)

plt.show()
```

---

## Internal Functions

The following functions are used internally but may be useful for advanced users:

### `_plot_point_forces(ax, point_forces, beam_length)`

Plots 3D arrows for point forces using lines and 3D cones. Point forces are rendered in red.

### `_plot_distributed_forces(ax, dist_forces, beam_length)`

Plots 3D arrows for distributed forces. Multiple arrows are drawn along the distributed load length, with all tails connected by a line. Distributed forces are rendered in green.

### `_plot_moments(ax, moments, beam_length)`

Plots 3D moment arrows as 3/4 arcs with cone tips. Moments are rendered in blue.

### `_plot_supports(ax, supports, beam_length)`

Plots supports as hollow circles (markers) with their 6-digit support type labels.

### `_create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=8)`

Creates a 3D cone mesh for arrow heads.

### `_create_moment_arc(center, moment_vec_beam, radius, arc_angle=3*π/2, num_points=30)`

Creates a 3/4 arc in a plane perpendicular to the moment vector.

### `_plot_stress_line(ax, loaded_beam, length, n_points=100)`

Plots the beam axis as a color-coded line based on von Mises stress using a Line3DCollection.

---

## Units

**Beamy is unit-agnostic.** Plotting functions use the same units as your input data. See the [Beam documentation](../beam.md#units) for details on unit consistency requirements.
