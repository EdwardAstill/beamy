# Stress Plotter

Reference documentation for plotting stress distributions on beam cross-sections at specific locations along the beam.

## Table of Contents

- [StressPlotter](#stressplotter)
  - [plot_stress_at](#plot_stress_at)

---

## StressPlotter

A class that plots stress distributions on the beam's cross-section at specific locations along the beam length. It integrates with `sectiony` to provide detailed stress visualization.

```python
from beamy import StressPlotter
# Or:
from beamy.analysis.stress_plotter import StressPlotter

class StressPlotter:
    """
    Plots stress distributions on the beam's cross-section at a specific location along the beam.
    """
    def __init__(self, loaded_beam: LoadedMember):
        """
        Initialize the stress plotter.
        
        Args:
            loaded_beam: LoadedMember object containing beam, loads, and analysis results
        """
```

### Initialization

```python
from beamy import LoadedMember, StressPlotter

# After creating a LoadedMember
lb = LoadedMember(beam, load_case)

# Create a StressPlotter
sp = StressPlotter(lb)
```

---

## plot_stress_at

Plots the stress distribution on the cross-section at a specific position along the beam.

```python
def plot_stress_at(
    self, 
    x_pos: float, 
    stress_type: StressType = "von_mises", 
    ax: Optional[plt.Axes] = None, 
    show: bool = True,
    cmap: str = "viridis",
    title: Optional[str] = None
) -> Optional[plt.Axes]:
    """
    Plot stress distribution on the section at a specific beam location x.
    
    Args:
        x_pos: Position along the beam length [0, L]
        stress_type: Type of stress to plot (e.g., 'von_mises', 'sigma_bending', etc.)
        ax: Matplotlib axes to plot on (creates new if None)
        show: Whether to show the plot immediately
        cmap: Colormap to use
        title: Optional title override
    
    Returns:
        The matplotlib axes object, or None if no geometry
    """
```

### Parameters

- `x_pos` (float): Position along the beam length where the stress distribution should be calculated. Should be in the range [0, L] where L is the beam length.
- `stress_type` (StressType): Type of stress to plot. Valid options:
  - `"von_mises"`: Von Mises equivalent stress (default)
  - `"sigma"`: Total normal stress (axial + bending)
  - `"sigma_axial"`: Normal stress due to axial force
  - `"sigma_bending"`: Normal stress due to bending moments
  - `"tau"`: Total shear stress magnitude
  - `"tau_shear"`: Shear stress due to transverse shear forces
  - `"tau_torsion"`: Shear stress due to torsion
- `ax` (matplotlib.axes.Axes, optional): Existing matplotlib axes to plot on. If None, creates a new figure and axes.
- `show` (bool): Whether to display the plot immediately using `plt.show()` (default: True)
- `cmap` (str): Colormap name to use for the stress visualization (default: "viridis")
- `title` (str, optional): Custom title for the plot. If None, uses a default title with the stress type and position.

### How It Works

1. **Force Extraction**: The method interpolates internal forces at the specified position `x_pos`:
   - Axial force (N) from `loaded_beam.axial()`
   - Shear forces (Vy, Vz) from `loaded_beam.shear("y")` and `loaded_beam.shear("z")`
   - Torsional moment (Mx) from `loaded_beam.torsion()`
   - Bending moments (My, Mz) from `loaded_beam.bending("z")` and `loaded_beam.bending("y")`
   
   **Note**: The bending moment mapping accounts for coordinate system differences:
   - `bending("y")` in beamy corresponds to moment about Z-axis (Mz) in sectiony
   - `bending("z")` in beamy corresponds to moment about Y-axis (My) in sectiony

2. **Stress Calculation**: Creates a `sectiony.Stress` object with the extracted forces and calculates the requested stress type.

3. **Visualization**: Uses `sectiony`'s plotting functionality to create a contour plot of the stress distribution on the cross-section.

### Example

```python
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedMember, StressPlotter
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
lb = LoadedMember(beam, lc)

# Find location of maximum von Mises stress
vm_results = lb.von_mises(points=200)
max_vm_idx = np.argmax(vm_results._values)
max_vm_x = vm_results._x[max_vm_idx]

print(f"Max Von Mises Stress at x = {max_vm_x:.2f}")

# Create stress plotter
sp = StressPlotter(lb)

# Plot von Mises stress at the critical location
sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="von_mises",
    title=f"Max Von Mises Stress Section (x={max_vm_x:.0f}mm)",
    cmap="plasma"
)

# Plot different stress types
sp.plot_stress_at(x_pos=2.5, stress_type="sigma_bending", cmap="RdBu_r")
sp.plot_stress_at(x_pos=2.5, stress_type="tau_shear", cmap="viridis")
```

### Plotting Multiple Sections

You can plot multiple sections on the same figure:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot at different locations
sp.plot_stress_at(x_pos=1.0, stress_type="von_mises", ax=axes[0], show=False, title="x=1.0m")
sp.plot_stress_at(x_pos=2.5, stress_type="von_mises", ax=axes[1], show=False, title="x=2.5m")
sp.plot_stress_at(x_pos=4.0, stress_type="von_mises", ax=axes[2], show=False, title="x=4.0m")

plt.tight_layout()
plt.show()
```

### Notes

- The section must have geometry defined (via `section.geometry`) for plotting to work
- Forces are interpolated from the analysis results, so accuracy depends on the number of points used in the analysis
- The method uses `np.interp` for linear interpolation of forces
- If `x_pos` is outside [0, L], a warning is printed but the plot will still attempt to proceed
- The plot includes a colorbar showing the stress scale
- Section outlines are drawn on top of the contour plot to hide any jagged edges from the grid-based masking

---

## Units

**Beamy is unit-agnostic.** The stress plotter uses the same units as your input data. The stress values will be in force/length² units (e.g., Pa if using N and m, or MPa if using kN and mm). See the [Beam documentation](../beam.md#units) for details on unit consistency requirements.
