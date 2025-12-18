from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List
import numpy as np
from ..core.support import validate_support_type

@dataclass
class Node:
    """A point in 3D space where members connect."""
    id: str
    position: np.ndarray
    support: Optional[str] = None
    connected_members: List[str] = field(default_factory=list)
    
    def __post_init__(self) -> None:
        if not isinstance(self.position, np.ndarray): self.position = np.array(self.position, dtype=float)
        if self.position.shape != (3,): raise ValueError(f"Node position must be 3D, got {self.position.shape}")
        if self.support: self.support = validate_support_type(self.support)
    
    def __repr__(self) -> str:
        pos = f"[{self.position[0]:.2f}, {self.position[1]:.2f}, {self.position[2]:.2f}]"
        return f"Node('{self.id}', pos={pos}{', support='+self.support if self.support else ''})"
