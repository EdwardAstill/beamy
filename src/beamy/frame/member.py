from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING
import numpy as np
from sectiony import Section
from ..core.material import Material

if TYPE_CHECKING:
    from .node import Node

@dataclass
class Member:
    """A beam element connecting two nodes in 3D space."""
    id: str
    start_node_id: str
    end_node_id: str
    section: Section
    material: Material
    orientation: np.ndarray
    releases: Optional[str] = None
    constraints: Optional[str] = None
    _start_node: Optional[Node] = None
    _end_node: Optional[Node] = None
    
    def __post_init__(self) -> None:
        if not isinstance(self.orientation, np.ndarray): self.orientation = np.array(self.orientation, dtype=float)
        if self.orientation.shape != (3,): raise ValueError(f"Orientation must be 3D vector, got {self.orientation.shape}")
        if np.linalg.norm(self.orientation) < 1e-10: raise ValueError(f"Member {self.id}: orientation cannot be zero")
        if self.releases and (len(self.releases) != 12 or not all(c in '01' for c in self.releases)):
            raise ValueError(f"Member {self.id}: releases must be 12-digit 0/1 string")
        if self.constraints and (len(self.constraints) != 12 or not all(c in "01" for c in self.constraints)):
            raise ValueError(f"Member {self.id}: constraints must be 12-digit 0/1 string")

    @property
    def length(self) -> float:
        if not (self._start_node and self._end_node): raise RuntimeError(f"Member {self.id}: nodes not set")
        return float(np.linalg.norm(self._end_node.position - self._start_node.position))
    
    @property
    def direction(self) -> np.ndarray:
        if not (self._start_node and self._end_node): raise RuntimeError(f"Member {self.id}: nodes not set")
        vec = self._end_node.position - self._start_node.position
        L = np.linalg.norm(vec)
        if L < 1e-10: raise ValueError(f"Member {self.id}: zero-length")
        return vec / L
    
    @property
    def local_y(self) -> np.ndarray:
        x = self.direction
        y_proj = self.orientation - np.dot(self.orientation, x) * x
        norm_y = np.linalg.norm(y_proj)
        if norm_y < 1e-10: raise ValueError(f"Member {self.id}: orientation is parallel to axis")
        return y_proj / norm_y
    
    @property
    def local_z(self) -> np.ndarray: return np.cross(self.direction, self.local_y)
    
    @property
    def transformation_matrix(self) -> np.ndarray: return np.array([self.direction, self.local_y, self.local_z])
    
    def set_nodes(self, start_node: Node, end_node: Node) -> None:
        self._start_node, self._end_node = start_node, end_node
    
    def __repr__(self) -> str: return f"Member('{self.id}', {self.start_node_id} -> {self.end_node_id})"
