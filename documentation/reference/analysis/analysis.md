# Analysis

Reference documentation for beam analysis functions and result classes.

## Table of Contents

- [Result](#result)
- [AnalysisResult](#analysisresult)
- [LoadedBeam](#loadedbeam)
- [solve_x_reactions](#solve_x_reactions)
- [solve_transverse_reactions](#solve_transverse_reactions)
- [get_all_loads](#get_all_loads)

---

## Result

Wrapper class for analysis results providing convenient access to (x, y) data pairs.

```python
from beamy import Result

class Result:
    def __init__(self, x: np.ndarray, values: np.ndarray):
        """
        Args:
            x: Array of x-coordinates
            values: Array of corresponding y-values
        """
```

### Properties

- `max` (float): Maximum value in the result
- `min` (float): Minimum value in the result
- `mean` (float): Mean value
- `range` (float): Range (max - min) of values

### Methods

#### `at(x_loc: float) -> float`

Interpolate the result at a specific x-location.

**Parameters:**
- `x_loc` (float): X-coordinate to evaluate at

**Returns:**
- `float`: Interpolated value at `x_loc`

### Usage

The `Result` class is iterable and can be indexed:

```python
# Iterate over (x, y) pairs
for x, y in result:
    print(f"x={x}, y={y}")

# Access by index
x, y = result[0]

# Get value at specific location
value = result.at(2.5)

# Get statistics
print(f"Max: {result.max}, Min: {result.min}, Mean: {result.mean}")
```

### Example

```python
import numpy as np
from beamy import Result

x = np.linspace(0, 5, 100)
y = np.sin(x)
result = Result(x, y)

print(result.max)      # Maximum value
print(result.at(2.5))  # Interpolated value at x=2.5
```

---

## AnalysisResult

Container for action, stress, and displacement results.

```python
from beamy import AnalysisResult

@dataclass
class AnalysisResult:
    _action: Result
    _stress: Result
    _displacement: Result
```

### Properties

- `action` (Result): Force or moment distribution (V, M, N, T)
- `stress` (Result): Stress distribution (sigma, tau)
- `displacement` (Result): Displacement distribution (w, u, theta)

### Example

```python
# Get analysis result from LoadedBeam
result = loaded_beam.shear(axis="z", points=100)

# Access action (shear force)
shear_force = result.action
print(shear_force.max)  # Maximum shear force

# Access stress (shear stress)
shear_stress = result.stress
print(shear_stress.at(2.5))  # Shear stress at x=2.5

# Access displacement
displacement = result.displacement
print(displacement.min)  # Minimum displacement
```

---

## LoadedBeam

Main analysis class that combines a beam with loads and provides analysis methods.

```python
from beamy import LoadedBeam

@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase
```

### Initialization

When a `LoadedBeam` is created, it automatically:
1. Solves for support reactions in all directions (x, y, z, rotations)
2. Combines applied loads with support reactions
3. Stores the complete load set for analysis

### Methods

#### `shear(axis: str, points: int = 100) -> AnalysisResult`

Calculate shear force and shear stress distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Shear force Fy and bending moment Mz
  - `"z"`: Shear force Fz and bending moment My
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (shear force), stress (shear stress), and displacement

#### `bending(axis: str, points: int = 100) -> AnalysisResult`

Calculate bending moment and bending stress distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Bending in x-z plane (moment My)
  - `"z"`: Bending in x-y plane (moment Mz)
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (bending moment), stress (bending stress), and displacement

#### `axial(points: int = 100) -> AnalysisResult`

Calculate axial force and axial stress distribution.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (axial force N), stress (axial stress σ = N/A), and displacement (axial displacement u)

#### `torsion(points: int = 100) -> AnalysisResult`

Calculate torsional moment and torsional stress distribution.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (torsional moment T), stress (torsional stress τ = T·r/J), and displacement (twist angle θ)

#### `deflection(axis: str, points: int = 100) -> Result`

Calculate transverse deflection distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Deflection in y-direction (v)
  - `"z"`: Deflection in z-direction (w)
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `Result`: Displacement distribution

#### `von_mises(points: int = 100) -> Result`

Calculate the Von Mises stress distribution along the beam.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `Result`: Von Mises stress distribution

**Note:** Uses a conservative superposition of maximum stress components:
- `sigma_vm = sqrt(sigma_max^2 + 3 * tau_max^2)`
- `sigma_max = |sigma_axial| + |sigma_bending_y| + |sigma_bending_z|`
- `tau_max = |tau_shear_y| + |tau_shear_z| + |tau_torsion|`

#### `plot(plot_stress: bool = False, plot_section: bool = True) -> None`

Plot the 3D beam diagram with loads and optional stress visualization.

**Parameters:**
- `plot_stress` (bool): If True, color the beam axis by von Mises stress (default: False)
- `plot_section` (bool): If True, draw the section outline at x=0 (default: True)

**Note:** This is a convenience method that calls `plot_beam_diagram()`. See the [Plotter documentation](plotter.md) for details.

### Example

```python
from beamy import Beam1D, Material, Section, Support, LoadCase, PointForce, LoadedBeam
import numpy as np

# Define beam
steel = Material(name="Steel", E=200e9, G=80e9)
section = Section(
    name="Rectangular",
    A=0.02, Iy=6.67e-5, Iz=1.67e-5, J=1e-6,
    y_max=0.05, z_max=0.1
)
beam = Beam1D(
    L=5.0,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),
        Support(x=5.0, type="111000")
    ]
)

# Define loads
lc = LoadCase(name="Point Load")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Create loaded beam and analyze
loaded_beam = LoadedBeam(beam, lc)

# Shear analysis
shear_result = loaded_beam.shear(axis="z", points=100)
print(f"Max shear: {shear_result.action.max}")  # force units
print(f"Max shear stress: {shear_result.stress.max}")  # force/length² units

# Bending analysis
bending_result = loaded_beam.bending(axis="z", points=100)
print(f"Max moment: {bending_result.action.max}")  # force×length units
print(f"Max bending stress: {bending_result.stress.max}")  # force/length² units

# Deflection
deflection = loaded_beam.deflection(axis="z", points=100)
print(f"Max deflection: {deflection.max}")  # length units
print(f"Deflection at midspan: {deflection.at(2.5)}")  # length units

# Axial analysis
axial_result = loaded_beam.axial(points=100)

# Torsion analysis
torsion_result = loaded_beam.torsion(points=100)

# Von Mises stress
von_mises = loaded_beam.von_mises(points=100)
print(f"Max von Mises stress: {von_mises.max}")  # force/length² units

# Plot beam diagram
loaded_beam.plot()
```

---

## solve_x_reactions

Solve for axial and torsional reactions using 1D finite element method.

```python
from beamy import solve_x_reactions

def solve_x_reactions(
    supports: List[Support],
    loads: LoadCase
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve for axial (Fx) and torsional (Mx) reactions.
    
    Args:
        supports: List of Support objects (must be sorted by x)
        loads: LoadCase containing applied loads
        
    Returns:
        Tuple of (d_x, d_rx) where:
        - d_x: Axial displacements at support nodes
        - d_rx: Torsional rotations at support nodes
    """
```

### Parameters

- `supports` (List[Support]): List of support conditions (should be sorted by x-position)
- `loads` (LoadCase): Load case containing applied loads

### Returns

- `Tuple[np.ndarray, np.ndarray]`: 
  - First array: Axial displacements (Ux) at each support node
  - Second array: Torsional rotations (Rx) at each support node

### Behavior

- Uses 1D finite element method with linear shape functions
- Automatically updates support reactions in the `Support.reactions` dictionary
- Reactions are stored with keys: `"Fx"` for axial reactions and `"Mx"` for torsional reactions

### Example

```python
from beamy import solve_x_reactions

# Supports should be sorted by x
supports_sorted = sorted(beam.supports, key=lambda s: s.x)

# Solve for reactions
d_x, d_rx = solve_x_reactions(supports_sorted, load_case)

# Reactions are automatically stored in support.reactions
for support in supports_sorted:
    print(f"At x={support.x}: Fx={support.reactions['Fx']}, Mx={support.reactions['Mx']}")
```

---

## solve_transverse_reactions

Solve for transverse reactions (shear and bending) using 1D finite element method with Hermite shape functions.

```python
from beamy import solve_transverse_reactions

def solve_transverse_reactions(
    beam: Beam1D,
    loads: LoadCase,
    axis: str = "z"
) -> np.ndarray:
    """
    Solve for transverse reactions (shear and bending moment).
    
    Args:
        beam: Beam1D object
        loads: LoadCase containing applied loads
        axis: Bending axis, either "y" or "z" (default: "z")
        
    Returns:
        Displacement vector [w1, theta1, w2, theta2, ...] for the specified axis
    """
```

### Parameters

- `beam` (Beam1D): Beam object with supports
- `loads` (LoadCase): Load case containing applied loads
- `axis` (str): Bending axis, either `"y"` or `"z"` (default: `"z"`)
  - `"y"`: Bending in x-y plane (shear Fy, moment Mz)
  - `"z"`: Bending in x-z plane (shear Fz, moment My)

### Returns

- `np.ndarray`: Displacement vector containing [w1, θ1, w2, θ2, ...] where:
  - `w`: Transverse displacement
  - `θ`: Rotation angle
  - Subscripts indicate support node indices

### Behavior

- Uses 1D finite element method with Hermite cubic shape functions
- Automatically updates support reactions in the `Support.reactions` dictionary
- For `axis="z"`: Reactions stored as `"Fz"` and `"My"`
- For `axis="y"`: Reactions stored as `"Fy"` and `"Mz"`

### Example

```python
from beamy import solve_transverse_reactions

# Solve for reactions in z-direction (vertical)
d_z = solve_transverse_reactions(beam, load_case, axis="z")

# Reactions are automatically stored
for support in beam.supports:
    print(f"At x={support.x}: Fz={support.reactions['Fz']}, My={support.reactions['My']}")

# Solve for reactions in y-direction (horizontal)
d_y = solve_transverse_reactions(beam, load_case, axis="y")
```

---

## get_all_loads

Get a sorted list of all loads (applied + reactions) for analysis.

```python
from beamy import get_all_loads

def get_all_loads(
    loads: LoadCase,
    beam: Beam1D
) -> List[Tuple[float, str, float]]:
    """
    Returns a sorted list of all loads and support reactions.
    
    Args:
        loads: LoadCase containing applied loads
        beam: Beam1D object with supports (reactions should be computed)
        
    Returns:
        List of (x, type, magnitude) tuples where:
        - x: Position along beam
        - type: Load type ("Fx", "Fy", "Fz", "Mx", "My", "Mz", "Rx", "Ry", "Rz", "RMx", "RMy", "RMz")
        - magnitude: Load magnitude
    """
```

### Parameters

- `loads` (LoadCase): Load case containing applied loads
- `beam` (Beam1D): Beam object with supports (reactions must be computed first)

### Returns

- `List[Tuple[float, str, float]]`: Sorted list of (x, type, magnitude) tuples

### Load Types

- **Applied Forces**: `"Fx"`, `"Fy"`, `"Fz"`
- **Applied Moments**: `"Mx"`, `"My"`, `"Mz"`
- **Reactions**: `"Rx"`, `"Ry"`, `"Rz"` (force reactions), `"RMx"`, `"RMy"`, `"RMz"` (moment reactions)

### Behavior

- Combines applied loads from the `LoadCase` with support reactions from `beam.supports`
- Loads at the same position and of the same type are summed
- Results are sorted by x-position

### Example

```python
from beamy import get_all_loads, LoadedBeam

# Create loaded beam (reactions are computed automatically)
loaded_beam = LoadedBeam(beam, load_case)

# Get all loads (applied + reactions)
all_loads = get_all_loads(load_case, beam)

# Print all loads
for x, load_type, magnitude in all_loads:
    print(f"At x={x:.2f}: {load_type} = {magnitude:.2f}")

# Filter for specific load type
vertical_forces = [(x, t, v) for x, t, v in all_loads if t in ("Fz", "Rz")]
```

---

## Analysis Workflow

Typical workflow for analyzing a beam:

1. **Define beam geometry and properties**
   ```python
   beam = Beam1D(L=5.0, material=steel, section=section, supports=[...])
   ```

2. **Define loads**
   ```python
   lc = LoadCase(name="My Loads")
   lc.add_point_force(...)
   lc.add_distributed_force(...)
   ```

3. **Create LoadedBeam** (automatically solves for reactions)
   ```python
   loaded_beam = LoadedBeam(beam, lc)
   ```

4. **Perform analysis**
   ```python
   shear_result = loaded_beam.shear(axis="z", points=100)
   bending_result = loaded_beam.bending(axis="z", points=100)
   deflection = loaded_beam.deflection(axis="z", points=100)
   ```

5. **Access results**
   ```python
   max_shear = shear_result.action.max
   max_stress = bending_result.stress.max
   max_deflection = deflection.max
   value_at_x = shear_result.action.at(2.5)
   ```

---

## Units

**Beamy is unit-agnostic.** All results use the same units as your inputs. See the [Beam documentation](../beam.md#units) for details on unit consistency requirements.

- **Action results** (force, moment): Same units as input forces/moments
- **Stress results**: Force / Length² (consistent with your input units)
- **Displacement results**: Length units (same as input length)
- **Rotation results**: Radians (dimensionless)

