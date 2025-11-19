# Beamy API Reference

Complete reference documentation for the `beamy` package.

## Table of Contents

- [Overview](#overview)
- [Beam Definitions](#beam-definitions)
  - [Support Strings](#support-strings)
  - [Material](#material)
  - [Section](#section)
  - [Beam1D](#beam1d)
- [Loads](#loads)
  - [PointLoad](#pointload)
  - [DistributedLoad](#distributedload)
  - [LoadCase](#loadcase)
- [Analysis](#analysis)
  - [BeamAnalysisResult](#beamanalysisresult)
  - [analyze_beam_simple_point_load](#analyze_beam_simple_point_load)
- [Utilities](#utilities)

---

## Overview

`beamy` is a lightweight 1D beam analysis package that provides:
- Beam definitions (material, section, supports)
- Point and distributed loads in local beam coordinates
- Closed-form static and Euler-Bernoulli bending analysis for a single beam

### Coordinate System

The beam uses a local coordinate system:
- **x**: Position along the beam axis (0 ≤ x ≤ L)
- **y**: Horizontal axis (perpendicular to beam)
- **z**: Vertical axis (perpendicular to beam)

Forces and moments follow this convention:
- **N**: Axial force (along x)
- **Vy, Vz**: Shear forces (along y, z)
- **My, Mz**: Bending moments (about y, z)
- **T**: Torsional moment (about x)

---

## Beam Definitions

### Support Strings

Beam supports are defined using 6-digit strings where each digit represents a constraint on a degree of freedom.

```python
from beamy import validate_support

support = "111000"  # 6-digit string
validate_support(support)  # Validates format
```

**Format:**
The support string is 6 digits representing constraints:
- **Position 1**: Ux (translation in x) - `0`=free, `1`=constrained
- **Position 2**: Uy (translation in y) - `0`=free, `1`=constrained
- **Position 3**: Uz (translation in z) - `0`=free, `1`=constrained
- **Position 4**: Rx (rotation about x) - `0`=free, `1`=constrained
- **Position 5**: Ry (rotation about y) - `0`=free, `1`=constrained
- **Position 6**: Rz (rotation about z) - `0`=free, `1`=constrained

**Common Support Types:**
- `"111111"`: Fully fixed (all DOFs constrained)
- `"111000"`: Pinned (translations constrained, rotations free)
- `"000000"`: Free (no constraints)
- `"110000"`: Roller in z-direction (Ux, Uy constrained; Uz and rotations free)
- `"100000"`: Roller in y and z-directions (only Ux constrained)

**Validation:**
The `validate_support()` function ensures the string is exactly 6 digits and contains only 0s and 1s.

**Example:**
```python
from beamy import validate_support

pinned = validate_support("111000")  # Valid
fixed = validate_support("111111")    # Valid
# validate_support("123")  # Raises ValueError
```

---

### Material

Material properties for the beam.

```python
from beamy import Material

@dataclass
class Material:
    name: str    # Material name
    E: float     # Young's modulus (Pa)
    G: float     # Shear modulus (Pa)
```

**Parameters:**
- `name` (str): Material identifier
- `E` (float): Young's modulus in Pascals (Pa)
- `G` (float): Shear modulus in Pascals (Pa)

**Example:**
```python
steel = Material(name="Steel", E=200e9, G=80e9)
```

---

### Section

Cross-section properties in local coordinates.

```python
from beamy import Section

@dataclass
class Section:
    name: str      # Section name
    A: float       # Area (m²)
    Iy: float      # Second moment of area about y-axis (m⁴)
    Iz: float      # Second moment of area about z-axis (m⁴)
    J: float       # Torsion constant (m⁴)
    y_max: float   # Distance to extreme fiber in +y direction (m)
    z_max: float   # Distance to extreme fiber in +z direction (m)
```

**Parameters:**
- `name` (str): Section identifier
- `A` (float): Cross-sectional area in square meters (m²)
- `Iy` (float): Second moment of area about local y-axis in m⁴
- `Iz` (float): Second moment of area about local z-axis in m⁴
- `J` (float): Torsion constant in m⁴
- `y_max` (float): Distance to extreme fiber in +y direction (m)
- `z_max` (float): Distance to extreme fiber in +z direction (m)

**Example:**
```python
rectangular = Section(
    name="Rectangular 100x200",
    A=0.02,           # 0.02 m²
    Iy=6.67e-5,        # m⁴
    Iz=1.67e-5,        # m⁴
    J=1e-6,            # m⁴
    y_max=0.05,        # 50 mm
    z_max=0.1          # 100 mm
)
```

---

### Beam1D

Straight prismatic beam along local x-axis.

```python
from beamy import Beam1D

@dataclass
class Beam1D:
    L: float           # Beam length (m)
    material: Material # Material properties
    section: Section   # Cross-section properties
    support_left: str   # Left support string (6 digits: Ux, Uy, Uz, Rx, Ry, Rz)
    support_right: str # Right support string (6 digits: Ux, Uy, Uz, Rx, Ry, Rz)
```

**Parameters:**
- `L` (float): Beam length in meters (m)
- `material` (Material): Material properties
- `section` (Section): Cross-section properties
- `support_left` (str): Support condition at x=0 as 6-digit string
- `support_right` (str): Support condition at x=L as 6-digit string

**Support String Format:**
Each support string must be exactly 6 digits (0s and 1s) representing constraints:
- Positions 1-3: Translations (Ux, Uy, Uz)
- Positions 4-6: Rotations (Rx, Ry, Rz)

The support strings are automatically validated when creating a `Beam1D` instance.

**Example:**
```python
beam = Beam1D(
    L=5.0,
    material=steel,
    section=rectangular,
    support_left="111000",   # Pinned: translations constrained, rotations free
    support_right="111000"   # Pinned: translations constrained, rotations free
)
```

---

## Loads

### PointLoad

Point load applied at a specific position along the beam.

```python
from beamy import PointLoad
import numpy as np

@dataclass
class PointLoad:
    x: float              # Position along beam (0 ≤ x ≤ L)
    force: np.ndarray     # Force vector [N, Vy, Vz]
    moment: np.ndarray    # Moment vector [T, My, Mz] (optional)
```

**Parameters:**
- `x` (float): Position along the beam axis (0 ≤ x ≤ L)
- `force` (np.ndarray): Force vector of shape (3,) containing [N, Vy, Vz]
  - `N`: Axial force (positive = tension)
  - `Vy`: Shear force in y-direction
  - `Vz`: Shear force in z-direction (positive = upward)
- `moment` (np.ndarray, optional): Moment vector of shape (3,) containing [T, My, Mz]
  - `T`: Torsional moment
  - `My`: Bending moment about y-axis
  - `Mz`: Bending moment about z-axis
  - Default: `np.zeros(3)`

**Example:**
```python
# 10 kN downward point load at midspan
load = PointLoad(
    x=2.5,
    force=np.array([0.0, 0.0, -10_000.0])  # Negative = downward
)
```

---

### DistributedLoad

Uniform distributed load per unit length over a segment.

```python
from beamy import DistributedLoad
import numpy as np

@dataclass
class DistributedLoad:
    x_start: float        # Start position (0 ≤ x_start ≤ L)
    x_end: float          # End position (x_start ≤ x_end ≤ L)
    w: np.ndarray         # Load intensity [n_x, v_y, v_z] per unit length
```

**Parameters:**
- `x_start` (float): Start position along beam (0 ≤ x_start ≤ L)
- `x_end` (float): End position along beam (x_start ≤ x_end ≤ L)
- `w` (np.ndarray): Load intensity vector of shape (3,) containing [n_x, v_y, v_z]
  - `n_x`: Axial load per unit length
  - `v_y`: Distributed load in y-direction per unit length
  - `v_z`: Distributed load in z-direction per unit length (positive = upward)

**Example:**
```python
# 5 kN/m downward distributed load over first half of beam
dist_load = DistributedLoad(
    x_start=0.0,
    x_end=2.5,
    w=np.array([0.0, 0.0, -5000.0])  # -5 kN/m downward
)
```

---

### LoadCase

Container for point and distributed loads.

```python
from beamy import LoadCase

@dataclass
class LoadCase:
    name: str                              # Load case name
    point_loads: List[PointLoad]         # List of point loads
    dist_loads: List[DistributedLoad]     # List of distributed loads
```

**Parameters:**
- `name` (str): Identifier for the load case
- `point_loads` (List[PointLoad]): List of point loads (default: empty list)
- `dist_loads` (List[DistributedLoad]): List of distributed loads (default: empty list)

**Methods:**

#### `add_point_load(load: PointLoad) -> None`
Add a point load to the load case.

**Parameters:**
- `load` (PointLoad): Point load to add

#### `add_distributed_load(load: DistributedLoad) -> None`
Add a distributed load to the load case.

**Parameters:**
- `load` (DistributedLoad): Distributed load to add

**Example:**
```python
lc = LoadCase(name="Dead Load + Live Load")
lc.add_point_load(PointLoad(x=2.5, force=np.array([0, 0, -10000])))
lc.add_distributed_load(DistributedLoad(
    x_start=0, x_end=5.0, w=np.array([0, 0, -2000])
))
```

---

## Analysis

### BeamAnalysisResult

Result object containing analysis outputs.

```python
from beamy import BeamAnalysisResult

@dataclass
class BeamAnalysisResult:
    beam: Beam1D
    reactions: Dict[str, Dict[str, float]]
    N: Callable[[float], float]    # Axial force function
    Vy: Callable[[float], float]   # Shear force in y
    Vz: Callable[[float], float]   # Shear force in z
    My: Callable[[float], float]  # Bending moment about y
    Mz: Callable[[float], float]  # Bending moment about z
    T: Callable[[float], float]   # Torsional moment
    w: Callable[[float], float]   # Deflection in z-direction
    theta: Callable[[float], float]  # Twist angle
```

**Attributes:**
- `beam` (Beam1D): The analyzed beam
- `reactions` (Dict[str, Dict[str, float]]): Support reactions
  - Format: `{"left": {"Vz": value, ...}, "right": {"Vz": value, ...}}`
- `N(x)` (Callable): Axial force at position x
- `Vy(x)` (Callable): Shear force in y-direction at position x
- `Vz(x)` (Callable): Shear force in z-direction at position x
- `My(x)` (Callable): Bending moment about y-axis at position x
- `Mz(x)` (Callable): Bending moment about z-axis at position x
- `T(x)` (Callable): Torsional moment at position x
- `w(x)` (Callable): Vertical deflection (displacement in z-direction) at position x
- `theta(x)` (Callable): Twist angle at position x

**Methods:**

#### `internal_forces(x: float) -> Dict[str, float]`
Get all internal forces at a given position.

**Parameters:**
- `x` (float): Position along beam

**Returns:**
- `Dict[str, float]`: Dictionary with keys "N", "Vy", "Vz", "My", "Mz", "T"

**Example:**
```python
forces = result.internal_forces(2.5)
# Returns: {"N": 0.0, "Vy": 0.0, "Vz": 5000.0, "My": 0.0, "Mz": 12500.0, "T": 0.0}
```

#### `normal_stress(x: float, y: float, z: float) -> float`
Calculate combined normal stress from axial force and bending moments.

**Parameters:**
- `x` (float): Position along beam
- `y` (float): Distance from neutral axis in y-direction (local coordinates)
- `z` (float): Distance from neutral axis in z-direction (local coordinates)

**Returns:**
- `float`: Combined normal stress (Pa) = σ_axial + σ_bending_y + σ_bending_z

**Formula:**
```
σ = N/A + My·z/Iy + Mz·y/Iz
```

**Example:**
```python
# Stress at midspan, top fiber
stress = result.normal_stress(x=2.5, y=0.0, z=0.1)
```

---

### analyze_beam_simple_point_load

Level 1 solver for simply supported beam with a single point load.

```python
from beamy import analyze_beam_simple_point_load

def analyze_beam_simple_point_load(
    beam: Beam1D,
    lc: LoadCase
) -> BeamAnalysisResult:
    ...
```

**Parameters:**
- `beam` (Beam1D): Beam to analyze (must be PINNED-PINNED)
- `lc` (LoadCase): Load case (must contain exactly one PointLoad)

**Returns:**
- `BeamAnalysisResult`: Analysis results

**Constraints:**
- Beam must have `support_left == "111000"` and `support_right == "111000"` (pinned supports)
- Load case must contain exactly one `PointLoad`
- Load case must not contain any `DistributedLoad` (not yet supported)
- Only vertical component Vz is considered (no axial, no torsion)

**Raises:**
- `ValueError`: If beam supports are not pinned-pinned
- `ValueError`: If load case doesn't contain exactly one point load
- `ValueError`: If load case contains distributed loads
- `ValueError`: If point load position is outside beam span

**Example:**
```python
result = analyze_beam_simple_point_load(beam, lc)

# Access reactions
print(result.reactions)
# {"left": {"Vz": 5000.0}, "right": {"Vz": 5000.0}}

# Query internal forces
shear = result.Vz(2.5)      # Shear at midspan
moment = result.Mz(2.5)     # Moment at midspan
deflection = result.w(2.5)   # Deflection at midspan

# Get all forces at a point
forces = result.internal_forces(2.5)
```

---

## Utilities

### heaviside

Heaviside step function.

```python
from beamy.utils import heaviside

def heaviside(x: float) -> float:
    """Simple Heaviside step (0 for x<0, 1 for x>=0)."""
```

**Parameters:**
- `x` (float): Input value

**Returns:**
- `float`: 0.0 if x < 0.0, else 1.0

**Example:**
```python
heaviside(-1.0)  # 0.0
heaviside(0.0)   # 1.0
heaviside(1.0)   # 1.0
```

---

### clip_to_span

Clamp a position value to the beam span [0, L].

```python
from beamy.utils import clip_to_span

def clip_to_span(x: float, L: float) -> float:
    """Clamp x to [0, L] for safety."""
```

**Parameters:**
- `x` (float): Position value
- `L` (float): Beam length

**Returns:**
- `float`: Clamped value in range [0, L]

**Example:**
```python
clip_to_span(-0.5, 5.0)  # 0.0
clip_to_span(3.0, 5.0)   # 3.0
clip_to_span(6.0, 5.0)   # 5.0
```

---

## Complete Example

```python
import numpy as np
from beamy import (
    Beam1D, Material, Section,
    LoadCase, PointLoad,
    analyze_beam_simple_point_load
)

# Define material
steel = Material(name="Steel", E=200e9, G=80e9)

# Define section (rectangular 100mm x 200mm)
section = Section(
    name="Rectangular",
    A=0.02,           # 0.02 m²
    Iy=6.67e-5,        # m⁴
    Iz=1.67e-5,        # m⁴
    J=1e-6,            # m⁴
    y_max=0.05,        # 50 mm
    z_max=0.1          # 100 mm
)

# Define beam (5 m long, simply supported)
beam = Beam1D(
    L=5.0,
    material=steel,
    section=section,
    support_left="111000",   # Pinned: translations constrained, rotations free
    support_right="111000"   # Pinned: translations constrained, rotations free
)

# Define load case
load = PointLoad(
    x=2.5,  # Midspan
    force=np.array([0.0, 0.0, -10_000.0])  # 10 kN downward
)

lc = LoadCase(name="Point Load at Midspan")
lc.add_point_load(load)

# Analyze
result = analyze_beam_simple_point_load(beam, lc)

# Print reactions
print("Reactions:", result.reactions)

# Print results at various points
for x in np.linspace(0.0, 5.0, 11):
    print(f"x={x:.2f} m: "
          f"Vz={result.Vz(x):.2f} N, "
          f"Mz={result.Mz(x):.2f} N·m, "
          f"w={result.w(x):.6e} m")
```

---

## Units

All physical quantities use SI base units:
- **Length**: meters (m)
- **Force**: Newtons (N)
- **Moment**: Newton-meters (N·m)
- **Stress**: Pascals (Pa)
- **Young's modulus**: Pascals (Pa)
- **Area**: square meters (m²)
- **Second moment of area**: meters to the fourth power (m⁴)

See `documentation/units.md` for more details.

