# Beam Definitions

Reference documentation for beam geometry, materials, sections, and supports.

## Table of Contents

- [Material](#material)
- [Section](#section)
- [Support](#support)
- [Beam1D](#beam1d)
- [validate_support_type](#validate_support_type)

---

## Material

Material properties for the beam.

**Note:** Beamy is unit-agnostic. Use consistent units throughout your analysis. See [Units](#units) for details.

```python
from beamy import Material

@dataclass
class Material:
    name: str    # Material name
    E: float     # Young's modulus (force/length²)
    G: float     # Shear modulus (force/length²)
    transparency: bool = False  # Used for plotting (affects alpha)
```

### Parameters

- `name` (str): Material identifier
- `E` (float): Young's modulus (force/length²). Must be consistent with your length and force units.
- `G` (float): Shear modulus (force/length²). Must be consistent with your length and force units.
- `transparency` (bool): Optional flag used for plotting visualization (default: False)

### Example

```python
# Example using SI units (Pa = N/m²)
steel = Material(name="Steel", E=200e9, G=80e9)

# Example using US units (psi = lbf/in²)
# steel = Material(name="Steel", E=29e6, G=11.5e6)
```

---

## Section

Cross-section properties in local beam coordinates.

**Note:** `Section` is imported from the `sectiony` library, not defined in `beamy`. The `sectiony` library provides section geometry and properties.

```python
from sectiony import Section
from sectiony.library import i_section, rectangular_section  # Common section types

# Section objects from sectiony have properties:
# - A: Area (length²)
# - Iy: Second moment of area about y-axis (length⁴)
# - Iz: Second moment of area about z-axis (length⁴)
# - J: Torsion constant (length⁴)
# - y_max: Distance to extreme fiber in +y direction (length)
# - z_max: Distance to extreme fiber in +z direction (length)
# - geometry: Geometry object containing shape definitions for plotting
```

### Required Properties

A `Section` object used with `Beam1D` must have the following properties:
- `A` (float): Cross-sectional area (length²)
- `Iy` (float): Second moment of area about local y-axis (length⁴)
- `Iz` (float): Second moment of area about local z-axis (length⁴)
- `J` (float): Torsion constant (length⁴)
- `y_max` (float): Distance to extreme fiber in +y direction (length)
- `z_max` (float): Distance to extreme fiber in +z direction (length)
- `geometry` (optional): Geometry object for visualization (used by plotting functions)

### Coordinate System

The beam uses a local coordinate system:
- **x**: Position along the beam axis (0 ≤ x ≤ L)
- **y**: Horizontal axis (perpendicular to beam)
- **z**: Vertical axis (perpendicular to beam)

### Example

```python
from sectiony.library import i_section

# I-beam section using sectiony library
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# Access properties (units depend on your input dimensions)
print(f"Area: {section.A}")
print(f"Iy: {section.Iy}")
print(f"Iz: {section.Iz}")
```

---

## Support

Support condition at a specific position along the beam.

```python
from beamy import Support

@dataclass
class Support:
    x: float
    type: str
    reactions: dict[str, float] = field(default_factory=lambda: {
        "Fx": 0.0, "Fy": 0.0, "Fz": 0.0,
        "Mx": 0.0, "My": 0.0, "Mz": 0.0,
    })
```

### Parameters

- `x` (float): Position along the beam axis where the support is located (0 ≤ x ≤ L)
- `type` (str): Support type as a 6-digit string (see [Support Strings](#support-strings))
- `reactions` (dict): Dictionary of support reactions (populated after analysis)

### Support Strings

The `type` parameter is a 6-digit string where each digit represents a constraint on a degree of freedom:

- **Position 1**: Ux (translation in x) - `0`=free, `1`=constrained
- **Position 2**: Uy (translation in y) - `0`=free, `1`=constrained
- **Position 3**: Uz (translation in z) - `0`=free, `1`=constrained
- **Position 4**: Rx (rotation about x) - `0`=free, `1`=constrained
- **Position 5**: Ry (rotation about y) - `0`=free, `1`=constrained
- **Position 6**: Rz (rotation about z) - `0`=free, `1`=constrained

### Common Support Types

- `"111111"`: Fully fixed (all DOFs constrained)
- `"111000"`: Pinned (translations constrained, rotations free)
- `"000000"`: Free (no constraints)
- `"110000"`: Roller in z-direction (Ux, Uy constrained; Uz and rotations free)
- `"100000"`: Roller in y and z-directions (only Ux constrained)

### Example

```python
# Pinned support at x=0
left_support = Support(x=0.0, type="111000")

# Fixed support at x=5.0
right_support = Support(x=5.0, type="111111")
```

---

## Beam1D

Straight prismatic beam along local x-axis.

```python
from beamy import Beam1D

@dataclass
class Beam1D:
    L: float           # Beam length (m)
    material: Material # Material properties
    section: Section   # Cross-section properties
    supports: List[Support]  # List of support conditions
```

### Parameters

- `L` (float): Beam length (length units)
- `material` (Material): Material properties
- `section` (Section): Cross-section properties
- `supports` (List[Support]): List of support conditions at various positions along the beam

### Validation

The beam automatically validates:
- Support strings are valid 6-digit strings
- Beam is supported in x, y, and z directions (at least one support constrains each translation)
- Beam is supported in rotation about x (at least one support constrains Rx)

### Example

```python
# Simply supported beam
beam = Beam1D(
    L=5.0,
    material=steel,
    section=rectangular,
    supports=[
        Support(x=0.0, type="111000"),   # Pinned at left end
        Support(x=5.0, type="111000")    # Pinned at right end
    ]
)

# Cantilever beam
cantilever = Beam1D(
    L=3.0,
    material=steel,
    section=rectangular,
    supports=[
        Support(x=0.0, type="111111")    # Fixed at left end
    ]
)
```

---

## validate_support_type

Validate a support type string format.

```python
from beamy import validate_support_type

def validate_support_type(support: str) -> str:
    """
    Validate node string format.
    
    Args:
        support: 6-digit string of 0s and 1s
        
    Returns:
        Validated node string
        
    Raises:
        ValueError: If node string is invalid
    """
```

### Parameters

- `support` (str): 6-digit string representing support constraints

### Returns

- `str`: The validated support string

### Raises

- `ValueError`: If the string is not exactly 6 digits or contains characters other than 0 and 1

### Example

```python
validate_support_type("111000")  # Returns "111000"
validate_support_type("111111")  # Returns "111111"
# validate_support_type("123")  # Raises ValueError
# validate_support_type("1110000")  # Raises ValueError
```

---

## Units

**Beamy is unit-agnostic.** You can use any consistent system of units (SI, US customary, etc.), as long as all quantities use the same base units throughout your analysis.

### Unit Consistency Requirements

All quantities must use consistent units:
- **Length**: Any length unit (m, ft, in, etc.)
- **Force**: Any force unit (N, lbf, kip, etc.)
- **Moment**: Force × Length (e.g., N·m, lbf·ft, kip·in)
- **Stress**: Force / Length² (e.g., Pa, psi, ksi)
- **Young's modulus**: Force / Length² (must match stress units)
- **Area**: Length²
- **Second moment of area**: Length⁴

### Example Unit Systems

**SI Units:**
- Length: meters (m)
- Force: Newtons (N)
- Stress/Modulus: Pascals (Pa = N/m²)
- Moment: Newton-meters (N·m)

**US Customary:**
- Length: inches (in) or feet (ft)
- Force: pounds-force (lbf) or kips (kip = 1000 lbf)
- Stress/Modulus: psi (lbf/in²) or ksi (kip/in²)
- Moment: lbf·in, lbf·ft, or kip·in, kip·ft

**Important:** All inputs for a single analysis must use the same unit system. Mixing units will produce incorrect results.

