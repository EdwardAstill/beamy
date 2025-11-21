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

```python
from beamy import Material

@dataclass
class Material:
    name: str    # Material name
    E: float     # Young's modulus (Pa)
    G: float     # Shear modulus (Pa)
```

### Parameters

- `name` (str): Material identifier
- `E` (float): Young's modulus in Pascals (Pa)
- `G` (float): Shear modulus in Pascals (Pa)

### Example

```python
steel = Material(name="Steel", E=200e9, G=80e9)
aluminum = Material(name="Aluminum", E=70e9, G=26e9)
```

---

## Section

Cross-section properties in local beam coordinates.

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

### Parameters

- `name` (str): Section identifier
- `A` (float): Cross-sectional area in square meters (m²)
- `Iy` (float): Second moment of area about local y-axis in m⁴
- `Iz` (float): Second moment of area about local z-axis in m⁴
- `J` (float): Torsion constant in m⁴
- `y_max` (float): Distance to extreme fiber in +y direction (m)
- `z_max` (float): Distance to extreme fiber in +z direction (m)

### Coordinate System

The beam uses a local coordinate system:
- **x**: Position along the beam axis (0 ≤ x ≤ L)
- **y**: Horizontal axis (perpendicular to beam)
- **z**: Vertical axis (perpendicular to beam)

### Example

```python
# Rectangular section: 100mm x 200mm
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

- `L` (float): Beam length in meters (m)
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

All physical quantities use SI base units:
- **Length**: meters (m)
- **Force**: Newtons (N)
- **Moment**: Newton-meters (N·m)
- **Stress**: Pascals (Pa)
- **Young's modulus**: Pascals (Pa)
- **Area**: square meters (m²)
- **Second moment of area**: meters to the fourth power (m⁴)

