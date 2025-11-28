# Beamy User Guide

Beamy is a lightweight Python package for 1D beam analysis, capable of handling static loads using Euler-Bernoulli beam theory. It supports axial, torsional, and transverse (shear/bending) analysis with 6 degrees of freedom per node.

**Note:** Beamy is unit-agnostic. Use consistent units throughout your analysis. See [Units](units.md) for details.

## Quick Start

Here is a complete example of setting up a simply supported beam with a point load.

```python
import numpy as np
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam
from sectiony.library import i_section

# 1. Define Properties
mat = Material(name="Steel", E=200e9, G=80e9)
sec = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)  # I-beam section

# 2. Create Beam
# Beam length (use consistent length units throughout)
L = 5.0
supports = [
    Support(x=0.0, type="111000"), # Pinned (Constrained x,y,z translations)
    Support(x=L, type="011000")    # Roller (Constrained y,z translations, free x)
]

beam = Beam1D(L=L, material=mat, section=sec, supports=supports)

# 3. Apply Loads
# Point force at mid-span (use consistent force units)
load = PointForce(
    point=np.array([2.5, 0.0, 0.0]),  # length units
    force=np.array([0.0, 0.0, -10000.0])  # force units (negative = downward)
)

case = LoadCase(name="Design Load")
case.add_point_force(load)

# 4. Solve
lb = LoadedBeam(beam, case)

# 5. Get Results
# Results use the same units as your inputs
print(f"Max Deflection: {lb.deflection('z').max:.6f}")  # length units
print(f"Max Bending Moment: {lb.bending('z').action.max:.2f}")  # force×length units
```

---

## Core Concepts

### 1. Material & Section
Define the physical properties of your beam.

```python
# Material: Young's Modulus (E) and Shear Modulus (G)
# Units must be force/length² (consistent with your length and force units)
steel = Material(name="Steel", E=200e9, G=77e9)  # Example: Pa if using m and N

# Section: Use sectiony library to create sections
from sectiony.library import i_section, rectangular_section

# I-beam section
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# Or rectangular section
# section = rectangular_section(width=0.1, height=0.2)
```

### 2. Supports
Supports are defined by a 6-digit string representing the 6 degrees of freedom (DOFs):
`[Ux, Uy, Uz, Rx, Ry, Rz]`
* `1` = Constrained (Fixed)
* `0` = Free

Common types:
* **Fixed/Clamped:** `"111111"`
* **Pinned:** `"111000"` (Translations fixed, rotations free)
* **Roller:** `"011000"` (Vertical/Lateral fixed, Axial free, Rotations free)

```python
from beamy import Support

# Support at x=0 fixed in all directions
s1 = Support(x=0.0, type="111111")
```

### 3. Loads
Loads are grouped into a `LoadCase`. You can add:
* `PointForce(point, force)`: Force vector `[Fx, Fy, Fz]` at `[x, y, z]`. Eccentric loads create moments.
* `Moment(x, moment)`: Moment vector `[Mx, My, Mz]` at location `x`.
* `DistributedForce(start_position, end_position, start_force, end_force)`: Linearly varying distributed force between two points.

```python
from beamy import LoadCase, PointForce, Moment, DistributedForce
import numpy as np

lc = LoadCase(name="Wind")
lc.add_point_force(PointForce(point=[2.0, 0, 0], force=[0, 100, 0]))
lc.add_moment(Moment(x=1.0, moment=[0, 0, 5000]))
lc.add_distributed_force(DistributedForce(
    start_position=np.array([0, 0, 0]),
    end_position=np.array([5, 0, 0]),
    start_force=np.array([0, -1000, 0]),
    end_force=np.array([0, -2000, 0])
))
```

### 4. Analysis
The `LoadedBeam` class automatically solves for reactions and internal forces upon initialization.

```python
lb = LoadedBeam(beam, load_case)
```

You can query specific behaviors using these methods:
* `lb.shear(axis, points=100)` - Shear force distribution
* `lb.bending(axis, points=100)` - Bending moment distribution
* `lb.axial(points=100)` - Axial force distribution
* `lb.torsion(points=100)` - Torsional moment distribution
* `lb.deflection(axis, points=100)` - Deflection distribution
* `lb.von_mises(points=100)` - Von Mises stress distribution along the beam

`axis` must be `"y"` or `"z"`.
* Bending about **z-axis** corresponds to loads in the **y-direction** (and vice versa depending on coordinate conventions, checking standard engineering axes is recommended). In `beamy`:
    * `bending("z")` -> Moment `Mz`, usually caused by loads in `y`.
    * `bending("y")` -> Moment `My`, usually caused by loads in `z`.
    *(Note: Check specific implementation details for sign conventions).*

### 5. Results
Analysis methods return an `AnalysisResult` object containing:
* `.action`: The internal force/moment (V, M, N, T).
* `.stress`: The calculated stress ($\sigma$ or $\tau$).
* `.displacement`: The relevant displacement ($w, u, \theta$).

Each of these properties is a `Result` object which wraps numpy arrays and provides helpers:

```python
res = lb.bending("z").action

res.max       # Maximum value
res.min       # Minimum value
res.at(1.5)   # Interpolated value at x=1.5m
res.mean      # Mean value

# Iterate over (x, value) pairs
for x, val in res:
    print(x, val)

# Von Mises stress (units: force/length²)
vm = lb.von_mises(points=100)
print(f"Max von Mises: {vm.max}")

# Plot the beam (convenience method)
lb.plot()  # Shows 3D diagram with loads
lb.plot(plot_stress=True)  # Shows 3D diagram with von Mises stress coloring
lb.plot(plot_stress=True, plot_section=False)  # Stress only, no section outline

# Or use plot_beam_diagram directly (more options)
from beamy.analysis.beam_plotter import plot_beam_diagram

plot_beam_diagram(lb, plot_stress=True, plot_section=True)
plot_beam_diagram(lb, plot_stress=True, save_path="beam_diagram.png")  # Save to file

# Plot section stress at specific locations
from beamy import StressPlotter

# Find location of maximum stress
vm_results = lb.von_mises(points=200)
max_vm_idx = np.argmax(vm_results._values)
max_vm_x = vm_results._x[max_vm_idx]

# Create stress plotter and plot stress at critical location
sp = StressPlotter(lb)
sp.plot_stress_at(
    x_pos=max_vm_x,
    stress_type="von_mises",
    title=f"Max Von Mises Stress (x={max_vm_x:.2f})",
    cmap="plasma"
)

# Plot different stress types at any location
sp.plot_stress_at(x_pos=2.5, stress_type="sigma_bending", cmap="RdBu_r")
sp.plot_stress_at(x_pos=2.5, stress_type="tau_shear", cmap="viridis")
```

---

## Visualization

Beamy provides two main plotting capabilities:

### 1. 3D Beam Diagram

The `plot_beam_diagram` function (or `lb.plot()` convenience method) creates a 3D visualization showing:
- Beam cross-section outline at x=0
- Beam axis (optionally colored by von Mises stress)
- Point forces as red arrows
- Distributed forces as green arrows with connecting lines
- Moments as blue arcs with cone tips
- Supports as hollow circles with labels

```python
from beamy.analysis.beam_plotter import plot_beam_diagram

# Basic plot
plot_beam_diagram(lb, plot_stress=False, plot_section=True)

# With stress coloring
plot_beam_diagram(lb, plot_stress=True, plot_section=True)

# Save to file instead of showing
plot_beam_diagram(lb, plot_stress=True, save_path="output.png")
```

### 2. Section Stress Plots

The `StressPlotter` class plots detailed stress distributions on the cross-section at specific locations:

```python
from beamy import StressPlotter

sp = StressPlotter(lb)

# Plot von Mises stress at a specific location
sp.plot_stress_at(x_pos=2.5, stress_type="von_mises")

# Available stress types:
# - "von_mises" (default)
# - "sigma" (total normal stress)
# - "sigma_axial" (axial stress)
# - "sigma_bending" (bending stress)
# - "tau" (total shear stress)
# - "tau_shear" (transverse shear)
# - "tau_torsion" (torsional shear)
```

For more details, see the [Beam Plotter](../reference/analysis/beam%20plotter.md) and [Section Plotter](../reference/analysis/section%20plotter.md) reference documentation.

