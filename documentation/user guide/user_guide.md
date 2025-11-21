# Beamy User Guide

Beamy is a lightweight Python package for 1D beam analysis, capable of handling static loads using Euler-Bernoulli beam theory. It supports axial, torsional, and transverse (shear/bending) analysis with 6 degrees of freedom per node.

## Quick Start

Here is a complete example of setting up a simply supported beam with a point load.

```python
import numpy as np
from beamy.beam import Beam1D, Material, Section, Support
from beamy.loads import LoadCase, PointForce
from beamy.analysis import LoadedBeam

# 1. Define Properties
mat = Material(name="Steel", E=200e9, G=80e9)
sec = Section(name="Rect", A=0.01, Iy=1e-6, Iz=1e-6, J=1e-6, y_max=0.1, z_max=0.1)

# 2. Create Beam
# A 5-meter beam supported at both ends
L = 5.0
supports = [
    Support(x=0.0, type="111000"), # Pinned (Constrained x,y,z translations)
    Support(x=L, type="011000")    # Roller (Constrained y,z translations, free x)
]

beam = Beam1D(L=L, material=mat, section=sec, supports=supports)

# 3. Apply Loads
# 10 kN downward force at mid-span
load = PointForce(
    point=np.array([2.5, 0.0, 0.0]), 
    force=np.array([0.0, 0.0, -10000.0])
)

case = LoadCase(name="Design Load")
case.add_point_force(load)

# 4. Solve
lb = LoadedBeam(beam, case)

# 5. Get Results
print(f"Max Deflection: {lb.deflection('z').max:.6f} m")
print(f"Max Bending Moment: {lb.bending('z').action.max:.2f} Nm")
```

---

## Core Concepts

### 1. Material & Section
Define the physical properties of your beam.

```python
# Material: Young's Modulus (E) and Shear Modulus (G) in Pascals
steel = Material(name="Steel", E=200e9, G=77e9)

# Section: Geometric properties
# A: Area (m^2)
# Iy, Iz: Second moments of area (m^4)
# J: Torsional constant (m^4)
# y_max, z_max: Max distances from centroid to fiber (for stress calc)
section = Section(name="IPE200", A=2.85e-3, Iy=1.94e-7, Iz=1.42e-5, J=6.98e-8, y_max=0.1, z_max=0.05)
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
from beamy.beam import Support

# Support at x=0 fixed in all directions
s1 = Support(x=0.0, type="111111")
```

### 3. Loads
Loads are grouped into a `LoadCase`. You can add:
* `PointForce(point, force)`: Force vector `[Fx, Fy, Fz]` at `[x, y, z]`. Eccentric loads create moments.
* `Moment(x, moment)`: Moment vector `[Mx, My, Mz]` at location `x`.
* `DistributedForce`: (Coming soon/In progress)

```python
lc = LoadCase(name="Wind")
lc.add_point_force(PointForce(point=[2.0, 0, 0], force=[0, 100, 0]))
```

### 4. Analysis
The `LoadedBeam` class automatically solves for reactions and internal forces upon initialization.

```python
lb = LoadedBeam(beam, load_case)
```

You can query specific behaviors using these methods:
* `lb.shear(axis)`
* `lb.bending(axis)`
* `lb.axial()`
* `lb.torsion()`
* `lb.deflection(axis)`
* `lb.von_mises()`

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
```

