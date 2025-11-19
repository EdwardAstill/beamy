# beam.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List

def validate_node(node: str) -> str:
    """
    Validate node string format.
    
    node string is 6 digits representing constraints:
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
        Validated node string
        
    Raises:
        ValueError: If node string is invalid
    """
    if not isinstance(node, str):
        raise ValueError(f"Support must be a string, got {type(node)}")
    if len(node) != 6:
        raise ValueError(f"node string must be exactly 6 digits, got '{node}' (length {len(node)})")
    if not all(c in '01' for c in node):
        raise ValueError(f"node string must contain only 0s and 1s, got '{node}'")
    return node

def validate_nodes(nodes: List[Node]):

    if not any(node.support[0] == '1' for node in nodes):
        raise ValueError("Beam must be supported in x direction at one or more nodes")
    if not any(node.support[1] == '1' for node in nodes):
        raise ValueError("Beam must be supported in y direction at one or more nodes")
    if not any(node.support[2] == '1' for node in nodes):
        raise ValueError("Beam must be supported in z direction at one or more nodes")
    if not any(node.support[3] == '1' for node in nodes):
        raise ValueError("Beam must be supported in rotation about x axis at one or more nodes")
    for node in nodes:
        validate_node(node.support)


@dataclass
class Material:
    name: str
    E: float  # Young's modulus (Pa)
    G: float  # Shear modulus (Pa)

@dataclass
class Section:
    """
    Cross-section properties in local coordinates.
    """
    name: str
    A: float   # area (m^2)
    Iy: float  # second moment about local y (m^4)
    Iz: float  # second moment about local z (m^4)
    J: float   # torsion constant (m^4)
    y_max: float  # distance to extreme fibre in +y (m)
    z_max: float  # distance to extreme fibre in +z (m)

@dataclass
class Node:
    """
    A node represents a point along the beam axis with support conditions.
    
    x: Position along beam axis (0 <= x <= L)
    support: 6-digit string representing constraints (validated using validate_support)
    """
    x: float
    support: str
    
    def __post_init__(self):
        """Validate the node string using validate_support."""
        self.support = validate_node(self.support)

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
    nodes: List[Node]

    def __post_init__(self):
        """Validate the nodes using validate_nodes."""
        validate_nodes(self.nodes)
        self.nodes = validate_nodes(self.nodes)
