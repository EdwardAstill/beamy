# Loads

Reference documentation for load definitions and load cases.

## Table of Contents

- [PointForce](#pointforce)
- [DistributedForce](#distributedforce)
- [Moment](#moment)
- [LoadCase](#loadcase)

---

## PointForce

Point force applied at a specific 3D location.

```python
from beamy import PointForce
import numpy as np

@dataclass
class PointForce:
    point: np.ndarray   # [x, y, z] - Position in beam coordinates
    force: np.ndarray   # [Fx, Fy, Fz] - Force vector
```

### Parameters

- `point` (np.ndarray): Position vector of shape (3,) containing [x, y, z]
  - `x`: Position along the beam axis (0 ≤ x ≤ L) (length units)
  - `y`: Horizontal offset from beam axis (length units)
  - `z`: Vertical offset from beam axis (length units)
- `force` (np.ndarray): Force vector of shape (3,) containing [Fx, Fy, Fz]
  - `Fx`: Axial force (positive = tension) (force units)
  - `Fy`: Shear force in y-direction (force units)
  - `Fz`: Shear force in z-direction (positive = upward) (force units)

### Coordinate System

Forces are defined in the local beam coordinate system:
- **x**: Along the beam axis
- **y**: Horizontal (perpendicular to beam)
- **z**: Vertical (perpendicular to beam)

### Example

```python
# Example using consistent units (e.g., meters and Newtons)
# Downward point load at midspan, on the beam axis
load1 = PointForce(
    point=np.array([2.5, 0.0, 0.0]),  # x=2.5 (length units)
    force=np.array([0.0, 0.0, -10_000.0])  # Fz=-10000 (force units, negative = downward)
)

# Horizontal force at x=1.0, offset in y-direction
load2 = PointForce(
    point=np.array([1.0, 0.1, 0.0]),  # x=1.0, y=0.1 (length units)
    force=np.array([0.0, 5_000.0, 0.0])  # Fy=5000 (force units)
)
```

---

## DistributedForce

Distributed force over a line segment in 3D space.

```python
from beamy import DistributedForce
import numpy as np

@dataclass
class DistributedForce:
    start_position: np.ndarray  # [x, y, z] - Start point
    end_position: np.ndarray    # [x, y, z] - End point
    start_force: np.ndarray    # [Fx, Fy, Fz] - Force per unit length at start
    end_force: np.ndarray       # [Fx, Fy, Fz] - Force per unit length at end
```

### Parameters

- `start_position` (np.ndarray): Start point vector of shape (3,) containing [x, y, z] (length units)
- `end_position` (np.ndarray): End point vector of shape (3,) containing [x, y, z] (length units)
- `start_force` (np.ndarray): Force intensity vector at start point, shape (3,) containing [Fx, Fy, Fz] per unit length (force/length)
- `end_force` (np.ndarray): Force intensity vector at end point, shape (3,) containing [Fx, Fy, Fz] per unit length (force/length)

### Behavior

The distributed force is linearly interpolated between the start and end points. The force intensity varies linearly along the segment.

### Example

```python
# Uniform distributed load (force per unit length)
# Example: 5 force_units/length_units downward
dist_load = DistributedForce(
    start_position=np.array([0.0, 0.0, 0.0]),  # length units
    end_position=np.array([2.5, 0.0, 0.0]),    # length units
    start_force=np.array([0.0, 0.0, -5_000.0]),  # force/length units
    end_force=np.array([0.0, 0.0, -5_000.0])     # force/length units (uniform)
)

# Linearly varying distributed load
varying_load = DistributedForce(
    start_position=np.array([0.0, 0.0, 0.0]),  # length units
    end_position=np.array([5.0, 0.0, 0.0]),    # length units
    start_force=np.array([0.0, 0.0, -10_000.0]),  # force/length units at start
    end_force=np.array([0.0, 0.0, 0.0])            # force/length units at end
)
```

---

## Moment

Applied moment at a specific position along the beam.

```python
from beamy import Moment
import numpy as np

@dataclass
class Moment:
    x: float
    moment: np.ndarray  # [T, My, Mz] - Moment vector
```

### Parameters

- `x` (float): Position along the beam axis (0 ≤ x ≤ L) (length units)
- `moment` (np.ndarray): Moment vector of shape (3,) containing [T, My, Mz] (force×length units)
  - `T`: Torsional moment (about x-axis)
  - `My`: Bending moment about y-axis
  - `Mz`: Bending moment about z-axis

### Example

```python
# Bending moment about z-axis at midspan (force×length units)
moment1 = Moment(
    x=2.5,  # length units
    moment=np.array([0.0, 0.0, 1000.0])  # force×length units
)

# Torsional moment (force×length units)
moment2 = Moment(
    x=1.0,  # length units
    moment=np.array([500.0, 0.0, 0.0])  # force×length units
)
```

---

## LoadCase

Container for point forces, distributed forces, and moments.

```python
from beamy import LoadCase

@dataclass
class LoadCase:
    name: str
    point_forces: List[PointForce] = field(default_factory=list)
    moments: List[Moment] = field(default_factory=list)
    dist_forces: List[DistributedForce] = field(default_factory=list)
    dist_force_resolution: int = 11  # Number of points to approximate distributed loads
```

### Parameters

- `name` (str): Identifier for the load case
- `point_forces` (List[PointForce]): List of point forces (default: empty list)
- `moments` (List[Moment]): List of applied moments (default: empty list)
- `dist_forces` (List[DistributedForce]): List of distributed forces (default: empty list)
- `dist_force_resolution` (int): Number of points used to discretize distributed loads (default: 11)

### Methods

#### `add_point_force(pf: PointForce) -> None`

Add a point force to the load case.

**Parameters:**
- `pf` (PointForce): Point force to add

#### `add_moment(m: Moment) -> None`

Add a moment to the load case.

**Parameters:**
- `m` (Moment): Moment to add

#### `add_distributed_force(df: DistributedForce) -> None`

Add a distributed force to the load case.

**Parameters:**
- `df` (DistributedForce): Distributed force to add

### Properties

The `LoadCase` class provides properties that automatically compute force and moment resultants:

#### `Fxs` -> List[Tuple[float, float]]

List of (x, Fx) pairs from all forces (explicit point forces + discretized distributed forces).

#### `Fys` -> List[Tuple[float, float]]

List of (x, Fy) pairs from all forces.

#### `Fzs` -> List[Tuple[float, float]]

List of (x, Fz) pairs from all forces.

#### `Mxs` -> List[Tuple[float, float]]

List of (x, Mx) pairs including:
- Torsion from eccentric forces (y*Fz - z*Fy)
- Explicit torsional moments

#### `Mys` -> List[Tuple[float, float]]

List of (x, My) pairs including:
- Bending from eccentric forces (z*Fx)
- Explicit bending moments about y-axis

#### `Mzs` -> List[Tuple[float, float]]

List of (x, Mz) pairs including:
- Bending from eccentric forces (-y*Fx)
- Explicit bending moments about z-axis

### Distributed Force Discretization

Distributed forces are automatically discretized into point forces when computing resultants. The `dist_force_resolution` parameter controls the number of discretization points (default: 11).

### Example

```python
lc = LoadCase(name="Dead Load + Live Load")

# Add point forces
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Add distributed load
lc.add_distributed_force(DistributedForce(
    start_position=np.array([0.0, 0.0, 0.0]),
    end_position=np.array([5.0, 0.0, 0.0]),
    start_force=np.array([0.0, 0.0, -2_000.0]),
    end_force=np.array([0.0, 0.0, -2_000.0])
))

# Add moment
lc.add_moment(Moment(
    x=1.0,
    moment=np.array([0.0, 0.0, 500.0])
))

# Access force resultants
print(lc.Fzs)  # List of (x, Fz) pairs
print(lc.Mzs)  # List of (x, Mz) pairs
```

---

## Units

**Beamy is unit-agnostic.** Use consistent units throughout your analysis. See the [Beam documentation](beam.md#units) for details on unit consistency requirements.

