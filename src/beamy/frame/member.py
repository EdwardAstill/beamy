from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING, Literal, List
import numpy as np
from sectiony import Section
from ..core.material import Material
from ..core.support import validate_support_type

if TYPE_CHECKING:
    from .node import Node


@dataclass
class Member:
    """A beam element defined by start/end positions in 3D space.
    
    The Frame automatically generates nodes and connectivity from member endpoints.
    Coincident endpoints are merged into single nodes.
    """
    id: str
    start: np.ndarray  # (3,) position vector
    end: np.ndarray    # (3,) position vector
    section: Section
    material: Material
    orientation: np.ndarray  # local Y direction
    element_type: Literal["beam", "truss", "cable"] = "beam"
    releases: Optional[str] = None
    constraints: Optional[str] = None
    _start_node: Optional[Node] = field(default=None, repr=False, compare=False)
    _end_node: Optional[Node] = field(default=None, repr=False, compare=False)
    
    def __post_init__(self) -> None:
        # Validate and convert positions
        if not isinstance(self.start, np.ndarray):
            object.__setattr__(self, 'start', np.array(self.start, dtype=float))
        if not isinstance(self.end, np.ndarray):
            object.__setattr__(self, 'end', np.array(self.end, dtype=float))
        if self.start.shape != (3,):
            raise ValueError(f"Member {self.id}: start must be 3D vector, got {self.start.shape}")
        if self.end.shape != (3,):
            raise ValueError(f"Member {self.id}: end must be 3D vector, got {self.end.shape}")
        
        # Check member length
        length = np.linalg.norm(self.end - self.start)
        if length < 1e-10:
            raise ValueError(f"Member {self.id}: zero-length member")
        
        # Validate orientation
        if not isinstance(self.orientation, np.ndarray):
            object.__setattr__(self, 'orientation', np.array(self.orientation, dtype=float))
        if self.orientation.shape != (3,):
            raise ValueError(f"Orientation must be 3D vector, got {self.orientation.shape}")
        if np.linalg.norm(self.orientation) < 1e-10:
            raise ValueError(f"Member {self.id}: orientation cannot be zero")
        
        # Validate element type
        if self.element_type not in ("beam", "truss", "cable"):
            raise ValueError(f"Member {self.id}: element_type must be 'beam', 'truss', or 'cable', got '{self.element_type}'")
        
        # Validate releases
        if self.releases and (len(self.releases) != 12 or not all(c in '01' for c in self.releases)):
            raise ValueError(f"Member {self.id}: releases must be 12-digit 0/1 string")
        
        # Validate constraints
        if self.constraints and (len(self.constraints) != 12 or not all(c in "01" for c in self.constraints)):
            raise ValueError(f"Member {self.id}: constraints must be 12-digit 0/1 string")

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.end - self.start))
    
    @property
    def direction(self) -> np.ndarray:
        vec = self.end - self.start
        L = np.linalg.norm(vec)
        if L < 1e-10:
            raise ValueError(f"Member {self.id}: zero-length")
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
    
    @property
    def start_node_id(self) -> str:
        """Node ID at the start of this member (set by Frame during construction)."""
        if self._start_node is None:
            raise RuntimeError(f"Member {self.id}: start node not set (Frame not built yet)")
        return self._start_node.id
    
    @property
    def end_node_id(self) -> str:
        """Node ID at the end of this member (set by Frame during construction)."""
        if self._end_node is None:
            raise RuntimeError(f"Member {self.id}: end node not set (Frame not built yet)")
        return self._end_node.id
    
    def __repr__(self) -> str:
        return f"Member('{self.id}', {tuple(self.start)} -> {tuple(self.end)})"
