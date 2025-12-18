# node.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List
import numpy as np

from ..setup.beam import validate_support_type


@dataclass
class Node:
    """
    A point in 3D space where members connect and where supports/loads can be applied.
    
    Attributes:
        id: Unique node identifier (e.g., "A", "N1", "base_left")
        position: 3D position vector [X, Y, Z] in global coordinates
        support: Optional support constraint string (6-digit, same format as Support.type)
        connected_members: List of member IDs connected to this node (auto-populated by Frame)
    """
    id: str
    position: np.ndarray
    support: Optional[str] = None
    connected_members: List[str] = field(default_factory=list)
    
    def __post_init__(self) -> None:
        """Validate the node after initialization."""
        # Validate position is a 3D array
        if not isinstance(self.position, np.ndarray):
            self.position = np.array(self.position, dtype=float)
        
        if self.position.shape != (3,):
            raise ValueError(f"Node position must be a 3D array [X, Y, Z], got shape {self.position.shape}")
        
        # Validate support string if provided
        if self.support is not None:
            self.support = validate_support_type(self.support)
    
    def __repr__(self) -> str:
        """String representation of the node."""
        pos_str = f"[{self.position[0]:.3f}, {self.position[1]:.3f}, {self.position[2]:.3f}]"
        support_str = f", support={self.support}" if self.support else ""
        return f"Node(id='{self.id}', position={pos_str}{support_str})"
