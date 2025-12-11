# beam.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

from sectiony import Section, Geometry

def validate_support_type(support: str) -> str:
    """
    Validate support string format.
    
    support string is 6 digits representing constraints:
    Position 1: Ux (translation in x) - 0=free, 1=constrained
    Position 2: Uy (translation in y) - 0=free, 1=constrained
    Position 3: Uz (translation in z) - 0=free, 1=constrained
    Position 4: Rx (rotation about x) - 0=free, 1=constrained
    Position 5: Ry (rotation about y) - 0=free, 1=constrained
    Position 6: Rz (rotation about z) - 0=free, 1=constrained
    
    Examples:
    - "111111": Fully fixed (all DOFs constrained)
    - "111000": Pinned (translations constrained, rotations free)
    - "000000": Free (no constraints)
    
    Args:
        support: 6-digit string of 0s and 1s
        
    Returns:
        Validated support string
        
    Raises:
        ValueError: If support string is invalid
    """
    if not isinstance(support, str):
        raise ValueError(f"Support must be a string, got {type(support)}")
    if len(support) != 6:
        raise ValueError(f"support string must be exactly 6 digits, got '{support}' (length {len(support)})")
    if not all(c in '01' for c in support):
        raise ValueError(f"support string must contain only 0s and 1s, got '{support}'")
    return support

def validate_support_pairs(supports: List[Support]):

    if not any(support.type[0] == '1' for support in supports):
        raise ValueError("Beam must be supported in x direction at one or more nodes")
    if not any(support.type[1] == '1' for support in supports):
        raise ValueError("Beam must be supported in y direction at one or more nodes")
    if not any(support.type[2] == '1' for support in supports):
        raise ValueError("Beam must be supported in z direction at one or more nodes")
    if not any(support.type[3] == '1' for support in supports):
        raise ValueError("Beam must be supported in rotation about x axis at one or more nodes")
    for support in supports:
        validate_support_type(support.type)


@dataclass
class Material:
    name: str
    E: float  # Young's modulus
    G: float  # Shear modulus
    Fy: Optional[float] = None  # Yield stress (needed for AISC checks)
    transparency: bool = False  # Used for plotting (affects alpha)


@dataclass
class Support:
    x: float
    type: str
    reactions: dict[str, float] = field(
        default_factory=lambda: {
            "Fx": 0.0, "Fy": 0.0, "Fz": 0.0,
            "Mx": 0.0, "My": 0.0, "Mz": 0.0,
        }
    )

    def __post_init__(self):
        self.type = validate_support_type(self.type)

@dataclass
class Beam1D:
    """
    Straight prismatic beam along local x in [0, L].
    
    Nodes are the points along the beam axis where the beam is supported.
    Each node contains a position (x) and support conditions (6-digit string).
    """
    L: float
    material: Material
    section: Section
    supports: List[Support]

    def __post_init__(self):
        """Validate the supports using validate_support_pairs."""
        validate_support_pairs(self.supports)
