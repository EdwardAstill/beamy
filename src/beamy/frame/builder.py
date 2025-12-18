"""
FrameBuilder: Coordinate-based frame definition.

Simplifies frame creation by allowing members to be defined by coordinates,
with nodes auto-generated from unique coordinate tuples.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
import numpy as np

from sectiony import Section
from ..core.material import Material
from .node import Node
from .member import Member
from .frame import Frame


# Type alias for coordinate tuple
Coord = Tuple[float, float, float]

# Tolerance for coordinate matching
COORD_TOL = 1e-9


def round_coord(coord: Coord, tol: float = COORD_TOL) -> Tuple[float, float, float]:
    """Round coordinate to tolerance for consistent hashing."""
    return (round(coord[0] / tol) * tol, round(coord[1] / tol) * tol, round(coord[2] / tol) * tol)


@dataclass
class MemberSpec:
    """Specification for a member before nodes are assigned."""
    id: str
    start: Coord
    end: Coord
    section: Section
    material: Material
    orientation: np.ndarray
    element_type: str = "beam"
    releases: Optional[str] = None
    constraints: Optional[str] = None


class FrameBuilder:
    """
    Builder for creating Frame objects from coordinate-based member definitions.
    
    Example:
        fb = FrameBuilder()
        fb.add("M1", (0, 0, 0), (1000, 0, 0), section, material)
        fb.add("M2", (1000, 0, 0), (1000, 0, 1000), section, material)
        fb.support_at((0, 0, 0), "111111")
        frame = fb.build()
    """
    
    def __init__(self, default_orientation: Tuple[float, float, float] = (0, 0, 1)):
        """
        Initialize the builder.
        
        Args:
            default_orientation: Default local Y-axis direction for members.
                                 (0, 0, 1) works for horizontal members.
        """
        self._members: List[MemberSpec] = []
        self._supports: Dict[Coord, str] = {}
        self._default_orientation = np.array(default_orientation, dtype=float)
    
    def add(
        self,
        member_id: str,
        start: Coord,
        end: Coord,
        section: Section,
        material: Material,
        orientation: Optional[Union[Tuple[float, float, float], np.ndarray]] = None,
        element_type: str = "beam",
        releases: Optional[str] = None,
        constraints: Optional[str] = None,
    ) -> "FrameBuilder":
        """
        Add a member defined by start and end coordinates.
        
        Args:
            member_id: Unique identifier for the member
            start: (x, y, z) start coordinate
            end: (x, y, z) end coordinate
            section: Cross-section properties
            material: Material properties
            orientation: Local Y-axis direction (defaults to builder's default)
            element_type: "beam", "truss", or "cable"
            releases: 12-digit string for DOF releases
            constraints: 12-digit string for DOF constraints
            
        Returns:
            self for method chaining
        """
        if orientation is None:
            orientation = self._default_orientation
        elif not isinstance(orientation, np.ndarray):
            orientation = np.array(orientation, dtype=float)
        
        self._members.append(MemberSpec(
            id=member_id,
            start=round_coord(start),
            end=round_coord(end),
            section=section,
            material=material,
            orientation=orientation,
            element_type=element_type,
            releases=releases,
            constraints=constraints,
        ))
        return self
    
    def support_at(self, coord: Coord, support_type: str) -> "FrameBuilder":
        """
        Add a support condition at a coordinate.
        
        Args:
            coord: (x, y, z) coordinate where support is applied
            support_type: 6-digit string (e.g., "111111" for fixed, "111000" for pinned)
            
        Returns:
            self for method chaining
        """
        self._supports[round_coord(coord)] = support_type
        return self
    
    def build(self) -> Frame:
        """
        Build the Frame from the defined members and supports.
        
        Nodes are auto-generated from unique coordinates with IDs N0, N1, N2, etc.
        
        Returns:
            Frame object ready for analysis
        """
        # Collect all unique coordinates
        coord_to_id: Dict[Coord, str] = {}
        node_counter = 0
        
        for spec in self._members:
            for coord in (spec.start, spec.end):
                if coord not in coord_to_id:
                    coord_to_id[coord] = f"N{node_counter}"
                    node_counter += 1
        
        # Create nodes
        nodes: List[Node] = []
        for coord, node_id in coord_to_id.items():
            support = self._supports[coord] if coord in self._supports else None
            nodes.append(Node(
                id=node_id,
                position=np.array(coord, dtype=float),
                support=support,
            ))
        
        # Create members
        members: List[Member] = []
        for spec in self._members:
            members.append(Member(
                id=spec.id,
                start_node_id=coord_to_id[spec.start],
                end_node_id=coord_to_id[spec.end],
                section=spec.section,
                material=spec.material,
                orientation=spec.orientation,
                element_type=spec.element_type,
                releases=spec.releases,
                constraints=spec.constraints,
            ))
        
        return Frame.from_nodes_and_members(nodes, members)
    
    def get_node_id_at(self, coord: Coord) -> Optional[str]:
        """
        Get the node ID that will be assigned to a coordinate (after build).
        
        This is useful for applying loads to specific nodes.
        Must be called after all members are added.
        
        Args:
            coord: (x, y, z) coordinate
            
        Returns:
            Node ID string or None if coordinate not found
        """
        rounded = round_coord(coord)
        # Rebuild the coord_to_id mapping
        coord_to_id: Dict[Coord, str] = {}
        node_counter = 0
        for spec in self._members:
            for c in (spec.start, spec.end):
                if c not in coord_to_id:
                    coord_to_id[c] = f"N{node_counter}"
                    node_counter += 1
        return coord_to_id[rounded] if rounded in coord_to_id else None
    
    def build_with_node_map(self) -> Tuple[Frame, Dict[Coord, str]]:
        """
        Build the Frame and return the coordinate-to-node-ID mapping.
        
        This is useful when you need to apply loads and need to know
        which node IDs correspond to which coordinates.
        
        Returns:
            (Frame, coord_to_node_id mapping)
        """
        # Build coord_to_id mapping
        coord_to_id: Dict[Coord, str] = {}
        node_counter = 0
        for spec in self._members:
            for coord in (spec.start, spec.end):
                if coord not in coord_to_id:
                    coord_to_id[coord] = f"N{node_counter}"
                    node_counter += 1
        
        # Create nodes
        nodes: List[Node] = []
        for coord, node_id in coord_to_id.items():
            support = self._supports[coord] if coord in self._supports else None
            nodes.append(Node(
                id=node_id,
                position=np.array(coord, dtype=float),
                support=support,
            ))
        
        # Create members
        members: List[Member] = []
        for spec in self._members:
            members.append(Member(
                id=spec.id,
                start_node_id=coord_to_id[spec.start],
                end_node_id=coord_to_id[spec.end],
                section=spec.section,
                material=spec.material,
                orientation=spec.orientation,
                element_type=spec.element_type,
                releases=spec.releases,
                constraints=spec.constraints,
            ))
        
        frame = Frame.from_nodes_and_members(nodes, members)
        return frame, coord_to_id
