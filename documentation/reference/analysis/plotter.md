# Plotter

Reference documentation for 3D visualization of beams and loads.

## Table of Contents

- [plot_extruded_section](#plot_extruded_section)
- [plot_loads](#plot_loads)

---

## plot_extruded_section

Plots a 3D extrusion of a 2D cross-section with optional load visualization.

```python
from beamy.plotter import plot_extruded_section

def plot_extruded_section(
    section_points: List[Tuple[float, float]],
    length: float,
    load_case: LoadCase = None
) -> None:
    """
    Plots a 3D extrusion of a 2D section.
    The section is defined in the YZ plane and extruded along the X axis.
    
    Args:
        section_points: List of (y, z) tuples defining the section vertices.
                       Must be a closed loop (first and last point same).
        length: Length of extrusion along the X axis.
        load_case: Optional LoadCase object containing forces to plot.
    """
```

### Parameters

- `section_points` (List[Tuple[float, float]]): List of (y, z) coordinate pairs defining the cross-section vertices. Must form a closed loop (first and last point should be the same).
- `length` (float): Length of the beam extrusion along the X axis (m)
- `load_case` (LoadCase, optional): Optional load case containing forces to visualize. If provided, point forces will be displayed as 3D arrows.

### Coordinate Mapping

The plot uses a specific coordinate mapping for visualization:
- **Plot X** ← Beam Z (Transverse Horizontal)
- **Plot Y** ← Beam X (Longitudinal)
- **Plot Z** ← Beam Y (Vertical)

This mapping provides a natural view where:
- The beam extends along the Y-axis (Plot Y = Beam X)
- The cross-section is viewed in the XZ plane (Plot X = Beam Z, Plot Z = Beam Y)

### Visual Features

- **Beam**: Silver-colored 3D mesh with semi-transparent faces
- **Point Forces**: Red 3D arrows with:
  - Arrow tip at the application point
  - Length proportional to force magnitude
  - Longest arrow is 1/5th of the beam length
  - 3D cone arrowhead for better visibility
- **Background**: White with light gray grid

### Example

```python
from beamy.plotter import plot_extruded_section
from beamy.loads import LoadCase, PointForce
import numpy as np

# Define triangular cross-section (y, z coordinates)
# Points form a closed loop: (0,0) -> (2,0) -> (1,1) -> (0,0)
points = ((0, 0), (2, 0), (1, 1), (0, 0))

# Create a load case
lc = LoadCase(name="Test Case")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.5, 0.2]),
    force=np.array([0, -1000, 500])
))
lc.add_point_force(PointForce(
    point=np.array([5.0, 1.0, 0.0]),
    force=np.array([-100, 0, 0])
))

# Plot the extruded section with loads
plot_extruded_section(points, length=5.0, load_case=lc)
```

### Notes

- The function creates a new matplotlib figure and displays it using `plt.show()`
- The view is set to elevation=20°, azimuth=-60° for a good default perspective
- Plot limits are automatically calculated with padding to ensure all geometry is visible
- At least 3 points are required to define a valid cross-section

---

## plot_beam_diagram

Plots a 3D beam diagram with section outline, loads, and optional stress visualization.

```python
from beamy.plotter import plot_beam_diagram

def plot_beam_diagram(
    loaded_beam: LoadedBeam,
    plot_stress: bool = False,
    plot_section: bool = True
) -> None:
    """
    Plots a 3D beam diagram with:
    - Optional 2D section outline on the YZ plane at x=0 (shear center at origin)
    - A line along the beam x-axis representing the beam length
    - Optional von Mises stress coloring on the beam axis
    - Point forces and moments as 3D arrows/arcs
    
    Args:
        loaded_beam: LoadedBeam object containing beam, loads, and analysis results
        plot_stress: If True, color the beam axis by von Mises stress
        plot_section: If True, draw the section outline at x=0
    """
```

### Parameters

- `loaded_beam` (LoadedBeam): LoadedBeam object containing beam, loads, and analysis results
- `plot_stress` (bool): If True, color the beam axis by von Mises stress (default: False)
- `plot_section` (bool): If True, draw the section outline at x=0 (default: True)

### Visual Features

- **Section Outline**: 2D cross-section drawn on the YZ plane at x=0, positioned with shear center at origin
- **Beam Axis**: Line along the beam length, optionally colored by von Mises stress
- **Point Forces**: Red 3D arrows showing force direction and magnitude
- **Moments**: Blue 3/4 arcs with cone tips showing moment direction
- **Stress Visualization**: When `plot_stress=True`, the beam axis is color-coded using a plasma colormap with a colorbar

### Coordinate Mapping

The plot uses a specific coordinate mapping for visualization:
- **Plot X** ← Beam Z (Transverse Horizontal)
- **Plot Y** ← Beam X (Longitudinal)
- **Plot Z** ← Beam Y (Vertical)

### Example

```python
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam
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
```

### Notes

- The function creates a new matplotlib figure and displays it using `plt.show()`
- The view is set to elevation=20°, azimuth=-60° for a good default perspective
- Plot limits are automatically calculated based on section size and beam length
- The section must have geometry defined (via `section.geometry`) for section plotting to work

---

## plot_loads

Plots loads on an existing 3D axes object.

```python
from beamy.plotter import plot_loads

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
- `beam_length` (float): Length of the beam (m), used to scale arrow lengths

### Behavior

- Plots point forces as 3D arrows
- Currently, moments and distributed forces are not visualized (functions are placeholders)

### Arrow Scaling

- Arrow lengths are proportional to force magnitude
- The longest arrow (corresponding to the maximum force magnitude) is set to `beam_length / 5`
- All other arrows scale proportionally

### Example

```python
import matplotlib.pyplot as plt
from beamy.plotter import plot_loads
from beamy.loads import LoadCase, PointForce
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

Plots 3D arrows for point forces using lines and 3D cones.

### `_create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=8)`

Creates a 3D cone mesh for arrow heads.

---

## Units

**Beamy is unit-agnostic.** Plotting functions use the same units as your input data. See the [Beam documentation](../beam.md#units) for details on unit consistency requirements.

