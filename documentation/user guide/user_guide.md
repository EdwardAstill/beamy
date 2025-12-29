# Beamy User Guide

Beamy is a lightweight Python package for 3D frame analysis with a convenient single-member wrapper (`LoadedMember`). It supports axial, torsional, and transverse (shear/bending) behavior with 6 degrees of freedom per node.

**Note:** Beamy is unit-agnostic. Use consistent units throughout your analysis. See [Units](units.md) for details.

Frames let you model connected members, joints, supports, trusses/cables, and loads across a structure. A single member is analyzed by building a 2-node / 1-member frame internally.

## Quick Start

Here is a complete example of setting up a simply supported member with a point load.

```python
import numpy as np
from sectiony.library import i_section
from beamy import Material, LoadCase, MemberPointForce, LoadedMember

# 1. Define Properties
mat = Material(name="Steel", E=200e9, G=80e9)
sec = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)  # I-beam section

# 2. Create member (use consistent length units throughout)
L = 5.0

# 3. Apply Loads
case = LoadCase(name="Design Load")
case.member_point_forces.append(
    MemberPointForce(
        member_id="M1",
        position=L / 2.0,
        force=np.array([0.0, 0.0, -10000.0]),
        coords="global",
        position_type="absolute",
    )
)

# 4. Solve
lb = LoadedMember(
    id="M1",
    start=np.array([0.0, 0.0, 0.0]),
    end=np.array([L, 0.0, 0.0]),
    section=sec,
    material=mat,
    orientation=np.array([0.0, 1.0, 0.0]),
    support_start="111100",
    support_end="011000",
    load_case=case,
)

# 5. Get Results
profile = lb.member_demand().actions(points=201)
print(f"Max bending moment |My|: {profile.bending_y.abs_max:.2f}")  # force×length units
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
The `LoadedMember` class automatically solves for reactions and internal forces upon initialization.

```python
lb = LoadedMember(beam, load_case)
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
plot_beam_diagram(lb, plot_stress=True, save_path="beam_diagram.svg")  # Save to file

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
plot_beam_diagram(lb, plot_stress=True, save_path="output.svg")
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

---

## Frames

Frames extend the same mechanics to assemblies of connected members in 3D. Each node has 6 DOFs `(Ux, Uy, Uz, Rx, Ry, Rz)`. Members connect at shared nodes and transfer forces/moments unless you add explicit end releases.

Key modules:
- Builder and core types in [src/beamy/frame/builder.py](src/beamy/frame/builder.py) and [src/beamy/frame/frame.py](src/beamy/frame/frame.py)
- Member definition in [src/beamy/frame/member.py](src/beamy/frame/member.py)
- Analysis and results in [src/beamy/frame/analysis.py](src/beamy/frame/analysis.py)
- Plotting in [src/beamy/viz/frame_plots.py](src/beamy/viz/frame_plots.py)

### Degrees of Freedom
- 6 DOF per node: translations `(Ux, Uy, Uz)` and rotations `(Rx, Ry, Rz)`.
- Supports constrain DOFs using 6-digit strings: `1` = fixed, `0` = free.
    - Fixed/clamped: `"111111"`
    - Pinned: `"111000"` (translations fixed, rotations free)
    - Roller examples: `"011000"`, `"101000"`, `"110000"` depending on the free translation.

### Nodes and Connectivity
- Members connect when they share the same coordinate (within builder rounding tolerance).
- If members cross without sharing a node at the intersection, they are not connected.
- The analysis pipeline auto-inserts real nodes where you attach loads/supports along members and splits members into segments for the solver. This happens inside `Frame.analyze(...)` (see solver expansion in [src/beamy/frame/analysis.py](src/beamy/frame/analysis.py)).

Best practice:
- Ensure any physical joint is represented by a shared node; split long members at connection points.
- Use member end releases to model pinned connections, not by leaving node rotations unconstrained at structural joints.

### Member End Releases and Constraints
Use releases to control moment transfer at member ends. These are 12-digit strings: first 6 digits apply to the start node DOFs, last 6 to the end node DOFs.

- Example: `"000111000111"` releases the three rotational DOFs at both ends (pin behavior for bending/torsion transfer).
- Constraints (also 12-digit strings) enforce end fixity beyond node supports when needed.

See `Member.releases` and `Member.constraints` in [src/beamy/frame/member.py](src/beamy/frame/member.py).

### Building Frames
Use `FrameBuilder` to define members by coordinates and add supports.

```python
import numpy as np
from sectiony.library import rhs
from beamy import Material
from beamy.frame import FrameBuilder
from beamy import LoadCase

steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
sec = rhs(b=0.150, h=0.150, t=0.006, r=0.0)

fb = FrameBuilder()

# Add two columns and a beam between them
fb.add("C1", (0.0, 0.0, 0.0), (0.0, 0.0, 3.0), sec, steel, orientation=(1,0,0))
fb.add("C2", (4.0, 0.0, 0.0), (4.0, 0.0, 3.0), sec, steel, orientation=(1,0,0))
fb.add("B1", (0.0, 0.0, 3.0), (4.0, 0.0, 3.0), sec, steel, orientation=(0,0,1))

# Base supports: fixed at both column bases
fb.support_at((0.0, 0.0, 0.0), "111111")
fb.support_at((4.0, 0.0, 0.0), "111111")

frame = fb.build()

# Loads: point load at mid-span of the beam (global coords)
loads = LoadCase("Service")
loads.add_member_point_load(
        member_id="B1",
        position=2.0,                      # absolute along member length (m)
        force=np.array([0.0, -10_000.0, 0.0]),
        moment=np.array([0.0, 0.0, 0.0]),
        coords="global",
        position_type="absolute",
)

frame.analyze(loads)
from beamy import plot_frame
plot_frame(frame, deformed=True)
```

### Supports
Add supports by coordinate via `FrameBuilder.support_at()`:

```python
fb.support_at((x, y, z), "111111")  # fixed
fb.support_at((x, y, z), "111000")  # pinned
```

At analysis time, those coordinates become nodes with the specified 6-DOF restraints. For truss/cable-only nodes (no beam elements attached), the solver auto-fixes rotations to avoid singular stiffness (these rotations have no physical stiffness contribution).

### Loads in Frames
Loads live in `LoadCase`:
- `add_nodal_force(node_id, force, coords)` and `add_nodal_moment(node_id, moment, coords)`
- `add_member_point_load(member_id, position, force, moment, coords, position_type)`
- `add_member_uniform_force(member_id, force_per_m, coords)` and distributed variants

Coordinate systems:
- `coords="global"`: vectors are in global XYZ.
- `coords="local"`: vectors are in the member’s local axes. Local axes are defined by the member direction and its `orientation` (local Y). See `Member.transformation_matrix` in [src/beamy/frame/member.py](src/beamy/frame/member.py).

Positions along members:
- `position_type="absolute"`: distance in the same units as your model.
- `position_type="relative"`: fraction of member length `[0,1]`.

During analysis setup, point loads/moments attached along members are converted to nodal loads at inserted split nodes. Distributed loads are split across segments. See [src/beamy/frame/analysis.py](src/beamy/frame/analysis.py).

### Analysis and Results
Call `frame.analyze(loads)` to run analysis. You can access:
- `frame.nodal_displacements`: per-node 6-DOF displacement vectors
- `frame.reactions`: per-supported node reaction vectors
- `frame.member_demand(member_id)`: solved member demand (`MemberActionProfile`) including axial, shear (`y/z`), torsion, and bending (`My/Mz`)

Example:
```python
demand = frame.member_demand("B1")
profile = demand.actions(points=801)
print("Max Mz:", profile.bending_z.abs_max)
print("Max Vy:", profile.shear_y.abs_max)
```

#### AISC Utilisation (Chapter F)
Use the plotting helper `plot_aisc_utilization(frame, ...)` (or run checks directly on `profile = frame.member_demand(member_id).actions(...)`).

### Visualization (Frames)
Plot functions accept a solved `Frame` (call `frame.analyze(load_case)` first for result-driven plots):
- `plot_frame(frame, ...)` – geometry (+ optional deformed overlay)
- `plot_deflection(frame, ...)` – colored deflections
- `plot_von_mises(frame, ...)` – stress coloring across members
- `plot_aisc_utilization(frame, ...)` – members colored by utilisation; one label per original member

Plot functions live in [src/beamy/viz/frame_plots.py](src/beamy/viz/frame_plots.py).

### Best Practices
- **Stability:** Constrain at least 6 independent global DOFs per load case (remove rigid-body modes). Don’t rely on incidental stiffness to stabilize the model.
- **Connectivity:** Split long members at real joints so intersecting members share nodes; otherwise they are not connected.
- **Pins vs Fixity:** Model pins using **member end releases** (e.g., release rotations at ends). Use node supports for boundary restraints, not for simulating member pins.
- **Truss/Cable-only nodes:** Rotations have no stiffness; allow the solver to auto-fix rotations at those nodes.
- **Coordinate discipline:** Be explicit about `coords="global"` vs `"local"` for loads; verify member local axes via `orientation`.
- **Sanity checks:** Verify reaction sums, symmetry (when expected), deformed shapes, and that utilisation makes sense.

---

## Frame Quick Start

```python
import numpy as np
from sectiony.library import rhs
from beamy import Material, LoadCase, plot_aisc_utilization
from beamy.frame import FrameBuilder

steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
sec = rhs(b=0.150, h=0.150, t=0.006, r=0.0)

fb = FrameBuilder()
fb.add("COL_L", (0,0,0), (0,0,2.5), sec, steel, orientation=(1,0,0))
fb.add("COL_R", (3,0,0), (3,0,2.5), sec, steel, orientation=(1,0,0))
fb.add("BEAM",  (0,0,2.5), (3,0,2.5), sec, steel, orientation=(0,0,1))

# Pin supports at bases
fb.support_at((0,0,0), "111000")
fb.support_at((3,0,0), "111000")

frame = fb.build()

loads = LoadCase("Demo")
loads.add_member_point_load(
        member_id="BEAM", position=1.5,
        force=np.array([0.0, -5_000.0, 0.0]),
        moment=np.array([0.0, 0.0, 0.0]),
        coords="global", position_type="absolute",
)

frame.analyze(loads)
plot_aisc_utilization(frame)
```

