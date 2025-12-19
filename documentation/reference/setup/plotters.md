````markdown
# Plotting Functions

Reference documentation for visualization of beams, frames, and structural results.

## Table of Contents

- [plot_beam_diagram](#plot_beam_diagram)
- [plot_frame](#plot_frame)
- [plot_deflection](#plot_deflection)
- [plot_von_mises](#plot_von_mises)
- [plot_results](#plot_results)
- [plot_member_diagrams](#plot_member_diagrams)
- [plot_analysis_results](#plot_analysis_results)
- [StressPlotter](#stressplotter)
- [plot_section](#plot_section)
- [plot_supports](#plot_supports)
- [plot_loads](#plot_loads)

---

## plot_beam_diagram

Plots a 3D beam diagram with section outline, loads, supports, and optional stress visualization.

```python
from beamy import plot_beam_diagram

def plot_beam_diagram(
    loaded_beam: LoadedMember,
    plot_stress: bool = False,
    plot_section: bool = True,
    save_path: Optional[str] = None
) -> None:
    """
    Plots a 3D beam diagram with:
    - Optional 2D section outline on the YZ plane at x=0
    - Beam axis line representing the beam length
    - Optional von Mises stress coloring on the beam axis
    - Point forces, distributed forces, and moments as 3D arrows/arcs
    - Supports as hollow circles with labels
    """
```

### Parameters

- `loaded_beam` (LoadedMember): LoadedMember object containing beam, loads, and analysis results
- `plot_stress` (bool): If True, color the beam axis by von Mises stress (default: False)
- `plot_section` (bool): If True, draw the section outline at x=0 (default: True)
- `save_path` (str, optional): If provided, saves the plot to this file path (high-quality .svg output)

### Visual Features

- **Section Outline**: 2D cross-section drawn on the YZ plane at x=0
- **Beam Axis**: Line along beam length, optionally colored by von Mises stress (plasma colormap)
- **Point Forces**: Red 3D arrows with proportional lengths
- **Distributed Forces**: Green 3D arrows along the load segment
- **Moments**: Blue 3/4 arcs with cone tips
- **Supports**: Hollow circles with 6-digit support type labels
- **Stress Visualization**: Colorbar for von Mises stress when plot_stress=True

### Coordinate Mapping

- Plot X ← Beam Z (transverse horizontal)
- Plot Y ← Beam X (longitudinal)
- Plot Z ← Beam Y (vertical)

### Example

```python
from beamy import Beam1D, Material, LoadedMember, LoadCase, PointForce
from beamy import plot_beam_diagram
from sectiony.library import i_section
import numpy as np

# Create beam
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)
beam = Beam1D(L=5.0, material=steel, section=section, supports=[
    Support(x=0.0, type="111000"),
    Support(x=5.0, type="111000")
])

# Create loads
lc = LoadCase(name="Point Load")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Plot
loaded_beam = LoadedMember(beam, lc)
plot_beam_diagram(loaded_beam, plot_stress=True, save_path="beam.svg")
```

---

## plot_frame

Plot the 3D frame geometry in wireframe style.

```python
from beamy import plot_frame

def plot_frame(
    loaded_frame: LoadedFrame,
    show_loads: bool = True,
    show_reactions: bool = True,
    show_member_ids: bool = True,
    show_node_ids: bool = True,
    deformed: bool = False,
    scale_factor: float = 1.0,
    save_path: Optional[str] = None
) -> None:
    """Plot frame geometry with optional deformation."""
```

### Parameters

- `loaded_frame` (LoadedFrame): Analyzed frame object
- `show_loads` (bool): Display applied load arrows (default: True)
- `show_reactions` (bool): Display reaction arrows (default: True)
- `show_member_ids` (bool): Label members (default: True)
- `show_node_ids` (bool): Label nodes (default: True)
- `deformed` (bool): Show deformed shape overlay (default: False)
- `scale_factor` (float): Deformation visualization scale (default: 1.0)
- `save_path` (str, optional): Save path for the plot

### Visual Features

- **Undeformed shape**: Solid black lines
- **Deformed shape**: Dashed blue lines (when deformed=True)
- **Supports**: Triangle markers at constrained nodes
- **Loads**: Red arrow quivers
- **Reactions**: Green arrow quivers

### Example

```python
loaded_frame = LoadedFrame(frame, loads)
plot_frame(loaded_frame, deformed=True, scale_factor=100, save_path="frame.svg")
```

---

## plot_deflection

Plot the deformed frame shape colored by displacement magnitude.

```python
from beamy import plot_deflection

def plot_deflection(
    loaded_frame: LoadedFrame,
    scale_factor: float = 1.0,
    points_per_member: int = 20,
    colormap: str = "viridis",
    show_undeformed: bool = True,
    show_colorbar: bool = True,
    save_path: Optional[str] = None
) -> None:
    """Plot deformed frame shape colored by displacement magnitude."""
```

### Parameters

- `loaded_frame` (LoadedFrame): Analyzed frame object
- `scale_factor` (float): Displacement scale multiplier (default: 1.0)
- `points_per_member` (int): Interpolation points per member (default: 20)
- `colormap` (str): Matplotlib colormap name (default: "viridis")
- `show_undeformed` (bool): Show original geometry (default: True)
- `show_colorbar` (bool): Display colorbar (default: True)
- `save_path` (str, optional): Save path

### Rendering

- Members discretized into segments
- Colors represent displacement magnitude: sqrt(Ux² + Uy² + Uz²)
- Undeformed geometry shown as faint dashed wireframe
- Colorbar shows displacement range

### Example

```python
plot_deflection(
    loaded_frame,
    scale_factor=100,
    colormap="plasma",
    save_path="deflection.svg"
)
```

---

## plot_von_mises

Plot the frame colored by Von Mises stress distribution.

```python
from beamy import plot_von_mises

def plot_von_mises(
    loaded_frame: LoadedFrame,
    points_per_member: int = 20,
    colormap: str = "jet",
    show_colorbar: bool = True,
    stress_limits: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> None:
    """Plot frame colored by Von Mises stress."""
```

### Parameters

- `loaded_frame` (LoadedFrame): Analyzed frame object
- `points_per_member` (int): Interpolation points per member (default: 20)
- `colormap` (str): Matplotlib colormap (default: "jet")
- `show_colorbar` (bool): Display colorbar (default: True)
- `stress_limits` (Tuple[float, float], optional): Fix color range (min, max)
- `save_path` (str, optional): Save path

### Stress Calculation

Von Mises stress: σ_vm = sqrt(σ² + 3τ²)

Where:
- σ = axial + bending stresses (conservative superposition)
- τ = shear + torsional stresses (conservative superposition)

### Color Interpretation

- **Warm colors (red/yellow)**: High stress
- **Cool colors (blue/green)**: Low stress

### Example

```python
plot_von_mises(
    loaded_frame,
    colormap="hot",
    stress_limits=(0, 345e6),  # 0 to yield stress
    save_path="stress.svg"
)
```

---

## plot_results

Unified plot showing deformed shape colored by analysis results.

```python
from beamy import plot_results

def plot_results(
    loaded_frame: LoadedFrame,
    result_type: str = "von_mises",
    deformed: bool = True,
    scale_factor: float = 1.0,
    points_per_member: int = 20,
    colormap: str = "jet",
    show_undeformed: bool = True,
    show_colorbar: bool = True,
    show_node_ids: bool = False,
    show_member_ids: bool = False,
    value_limits: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> None:
    """Plot frame colored by specified result type."""
```

### Parameters

- `loaded_frame` (LoadedFrame): Analyzed frame object
- `result_type` (str): Type of result to display (default: "von_mises")
  - "von_mises": Von Mises stress
  - "deflection": Displacement magnitude
  - "axial_stress": Axial stress
  - "bending_stress": Combined bending stress
  - "shear_stress": Combined shear stress
  - "axial_force": Axial force N
  - "shear_force": Combined shear force
  - "bending_moment": Combined bending moment
  - "torsion": Torsional moment
- `deformed` (bool): Plot on deformed geometry (default: True)
- `scale_factor` (float): Deformation scale (default: 1.0)
- `points_per_member` (int): Interpolation points (default: 20)
- `colormap` (str): Matplotlib colormap (default: "jet")
- `show_undeformed` (bool): Show original geometry (default: True)
- `show_colorbar` (bool): Display colorbar (default: True)
- `show_node_ids` (bool): Label nodes (default: False)
- `show_member_ids` (bool): Label members (default: False)
- `value_limits` (Tuple[float, float], optional): Color range limits
- `save_path` (str, optional): Save path

### Example

```python
# Von Mises stress on deformed shape
plot_results(
    loaded_frame,
    result_type="von_mises",
    deformed=True,
    scale_factor=50,
    save_path="results.svg"
)

# Bending moment distribution (use diverging colormap for +/- values)
plot_results(
    loaded_frame,
    result_type="bending_moment",
    deformed=False,
    colormap="RdBu_r",
    save_path="moment.svg"
)
```

---

## plot_member_diagrams

Plot internal force diagrams (N, V, M, T) for a specific member.

```python
from beamy import plot_member_diagrams

def plot_member_diagrams(
    loaded_frame: LoadedFrame,
    member_id: str,
    save_path: Optional[str] = None
) -> None:
    """Plot 2D internal force diagrams for a member."""
```

### Parameters

- `loaded_frame` (LoadedFrame): Analyzed frame object
- `member_id` (str): Member ID to plot
- `save_path` (str, optional): Save path

### Output

Creates a 2×2 subplot figure showing:

- **Top-left**: Axial force N(x)
- **Top-right**: Shear forces Vy(x) and Vz(x)
- **Bottom-left**: Bending moments My(x) and Mz(x)
- **Bottom-right**: Torsion T(x)

All diagrams use the member's local x-axis (0 to L).

### Example

```python
plot_member_diagrams(loaded_frame, "beam1", save_path="diagrams.svg")
```

---

## plot_analysis_results

Plots analysis results (shear, moment, deflection, axial/torsion) as 2D line graphs.

```python
from beamy import plot_analysis_results

def plot_analysis_results(
    loaded_beam: LoadedMember,
    save_path: Optional[str] = None,
    show: bool = True,
    points: int = 100
) -> None:
    """Plot analysis results as 2D line graphs."""
```

### Parameters

- `loaded_beam` (LoadedMember): Analyzed beam object
- `save_path` (str, optional): Save path
- `show` (bool): Display the plot (default: True)
- `points` (int): Number of sample points (default: 100)

### Output

Creates multi-panel plots showing:
- Shear force distributions
- Bending moment distributions
- Deflection distributions
- Axial and torsional distributions

### Example

```python
plot_analysis_results(loaded_beam, save_path="analysis.svg")
```

---

## StressPlotter

Utility class for advanced stress visualization.

```python
from beamy import StressPlotter

class StressPlotter:
    """Helper class for plotting stress distributions on frame members."""
```

### Methods

#### `plot_axial_stress(loaded_frame: LoadedFrame, colormap: str = "RdBu_r", ...) -> None`

Plot axial stress distribution on frame.

#### `plot_bending_stress(loaded_frame: LoadedFrame, component: str = "z", ...) -> None`

Plot bending stress distribution (from moment).

**Parameters:**
- `component` (str): Bending component, either "y" or "z" (default: "z")

#### `plot_combined_stress(loaded_frame: LoadedFrame, method: str = "von_mises", ...) -> None`

Plot combined stress using specified method.

**Parameters:**
- `method` (str): Combination method, options: "von_mises", "max", "tresca"

### Example

```python
from beamy import StressPlotter

plotter = StressPlotter()
plotter.plot_combined_stress(loaded_frame, method="von_mises", save_path="stress.svg")
```

---

## plot_section

Plot the cross-section geometry.

```python
from beamy import plot_section

def plot_section(
    section: Section,
    ax: Optional[matplotlib.axes.Axes] = None,
    show: bool = True
) -> Optional[matplotlib.axes.Axes]:
    """Plot the cross-section geometry."""
```

### Parameters

- `section` (Section): Section object to plot (from sectiony library)
- `ax` (matplotlib.axes.Axes, optional): Axes to plot on
- `show` (bool): Call plt.show() (default: True)

### Returns

- `matplotlib.axes.Axes`: The axes object

### Example

```python
from sectiony.library import i_section
from beamy import plot_section

section = i_section(d=0.3, b=0.2, tf=0.015, tw=0.010, r=0.0)
plot_section(section)
```

---

## plot_supports

Plots the beam supports in 2D.

```python
from beamy import plot_supports

def plot_supports(
    supports: List[Support],
    beam_length: float,
    unit: str = "m",
    save_path: Optional[str] = None,
    show: bool = True
) -> None:
    """Plot beam supports as a 2D diagram."""
```

### Parameters

- `supports` (List[Support]): List of support objects
- `beam_length` (float): Beam length
- `unit` (str): Unit label for position (default: "m")
- `save_path` (str, optional): Save path
- `show` (bool): Display the plot (default: True)

### Visual Features

- Beam shown as a straight line
- Supports marked as dots and labeled with support type

### Example

```python
from beamy import Support, plot_supports

supports = [
    Support(x=0.0, type="111000"),
    Support(x=5.0, type="111000")
]
plot_supports(supports, beam_length=5.0)
```

---

## plot_loads

Plots loads on a 3D axes object.

```python
from beamy import plot_loads

def plot_loads(
    ax: matplotlib.axes.Axes3D,
    load_case: LoadCase,
    beam_length: float
) -> None:
    """Plot loads on a 3D axes object."""
```

### Parameters

- `ax` (matplotlib.axes.Axes3D): 3D matplotlib axes
- `load_case` (LoadCase): Load case containing forces
- `beam_length` (float): Beam length for scaling

### Rendering

- **Point forces**: Red 3D arrows
- **Distributed forces**: Green 3D arrows with connecting line
- **Moments**: Blue 3/4 arcs with cone tips
- Arrow lengths proportional to magnitude

### Arrow Scaling

Maximum arrow length = beam_length / 5, other arrows scale proportionally.

### Example

```python
import matplotlib.pyplot as plt
from beamy import plot_loads, LoadCase, PointForce
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

lc = LoadCase(name="Forces")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

plot_loads(ax, lc, beam_length=5.0)
plt.show()
```

---

## Units

**Beamy is unit-agnostic.** Plotting functions use the same units as your input data. Ensure consistent units throughout your analysis.

### Common Unit Systems

**SI Units:**
- Length: meters (m)
- Force: Newtons (N)
- Stress: Pascals (Pa)

**US Customary:**
- Length: inches (in)
- Force: pounds-force (lbf) or kips (kip)
- Stress: psi or ksi

---

## Plotting Best Practices

### Color Map Selection

- **Jet**: Good for most analyses, full color range
- **Viridis**: Perceptually uniform, colorblind-friendly
- **Hot**: Good for stress (warm = high)
- **RdBu_r**: Diverging colormap for +/- moment diagrams
- **Plasma**: Good for deflection visualization

### Scale Factors

For stress visualization:
- Use `stress_limits` to fix colorbar at material yield stress for consistent interpretation
- Use `show_undeformed=True` to see original geometry for context

For deflection visualization:
- Use `scale_factor` to exaggerate small deflections for visibility
- Typical range: 50-500× depending on structure size and deflection

### File Formats

- **.svg**: Recommended for technical documents (vector format, scalable)
- **.png**: Good for presentations (raster format)
- **.pdf**: Good for publications

All plotting functions support `save_path` parameter for exporting high-quality images.

````
