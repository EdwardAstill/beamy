# member.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING
import numpy as np

from sectiony import Section
from ..setup.beam import Material

if TYPE_CHECKING:
    from .node import Node


@dataclass
class Member:
    """
    A beam element connecting two nodes in 3D space.
    
    Attributes:
        id: Unique member identifier (e.g., "M1", "column_1", "beam_AB")
        start_node_id: ID of the node at the start of the member
        end_node_id: ID of the node at the end of the member
        section: Cross-section properties (from sectiony library)
        material: Material properties (E, G, Fy)
        orientation: A 3D vector in global coordinates defining the local y-axis direction
        releases: Optional member end releases as 12-digit string [start 6 DOFs, end 6 DOFs]
    """
    id: str
    start_node_id: str
    end_node_id: str
    section: Section
    material: Material
    orientation: np.ndarray
    releases: Optional[str] = None
    
    # Private cache for computed properties (set by Frame)
    _start_node: Optional[Node] = None
    _end_node: Optional[Node] = None
    
    def __post_init__(self) -> None:
        """Validate the member after initialization."""
        # Validate orientation is a 3D array
        if not isinstance(self.orientation, np.ndarray):
            self.orientation = np.array(self.orientation, dtype=float)
        
        if self.orientation.shape != (3,):
            raise ValueError(f"Member orientation must be a 3D array, got shape {self.orientation.shape}")
        
        # Validate orientation is not zero vector
        if np.linalg.norm(self.orientation) < 1e-10:
            raise ValueError(f"Member {self.id}: orientation vector cannot be zero")
        
        # Validate releases if provided
        if self.releases is not None:
            if not isinstance(self.releases, str):
                raise ValueError(f"Member {self.id}: releases must be a string")
            if len(self.releases) != 12:
                raise ValueError(f"Member {self.id}: releases must be 12 digits, got {len(self.releases)}")
            if not all(c in '01' for c in self.releases):
                raise ValueError(f"Member {self.id}: releases must contain only 0s and 1s")
    
    @property
    def length(self) -> float:
        """Member length computed from node positions."""
        if self._start_node is None or self._end_node is None:
            raise RuntimeError(f"Member {self.id}: nodes not set. Member must be part of a Frame.")
        return float(np.linalg.norm(self._end_node.position - self._start_node.position))
    
    @property
    def direction(self) -> np.ndarray:
        """Unit vector from start to end node (local x-axis)."""
        if self._start_node is None or self._end_node is None:
            raise RuntimeError(f"Member {self.id}: nodes not set. Member must be part of a Frame.")
        vec = self._end_node.position - self._start_node.position
        length = np.linalg.norm(vec)
        if length < 1e-10:
            raise ValueError(f"Member {self.id}: zero-length member (start and end nodes at same position)")
        return vec / length
    
    @property
    def local_y(self) -> np.ndarray:
        """Unit vector for local y-axis (from orientation, orthogonalized to x)."""
        x_local = self.direction
        
        # Project orientation vector onto plane perpendicular to x_local (Gram-Schmidt)
        y_proj = self.orientation - np.dot(self.orientation, x_local) * x_local
        
        # Check if orientation is parallel to member axis
        norm_y = np.linalg.norm(y_proj)
        if norm_y < 1e-10:
            raise ValueError(
                f"Member {self.id}: orientation vector is parallel to member axis. "
                "Choose a different orientation vector."
            )
        
        return y_proj / norm_y
    
    @property
    def local_z(self) -> np.ndarray:
        """Unit vector for local z-axis (x × y, right-hand rule)."""
        return np.cross(self.direction, self.local_y)
    
    @property
    def transformation_matrix(self) -> np.ndarray:
        """
        3×3 rotation matrix from local to global coordinates.
        
        Each row contains a local axis expressed in global coordinates:
        Row 0: local x-axis (direction)
        Row 1: local y-axis
        Row 2: local z-axis
        """
        T = np.array([
            self.direction,
            self.local_y,
            self.local_z
        ])
        return T
    
    def set_nodes(self, start_node: Node, end_node: Node) -> None:
        """
        Set the node references for this member.
        Called by Frame after construction.
        """
        self._start_node = start_node
        self._end_node = end_node
    
    def __repr__(self) -> str:
        """String representation of the member."""
        return f"Member(id='{self.id}', {self.start_node_id} → {self.end_node_id})"
