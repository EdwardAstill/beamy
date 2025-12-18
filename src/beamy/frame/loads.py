# loads.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List
import numpy as np


@dataclass
class NodalForce:
    """
    Force applied directly at a node (in global coordinates).
    
    Attributes:
        node_id: Target node ID
        force: Force vector [FX, FY, FZ] in global coordinates
    """
    node_id: str
    force: np.ndarray
    
    def __post_init__(self) -> None:
        """Validate the nodal force."""
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)
        if self.force.shape != (3,):
            raise ValueError(f"NodalForce must be a 3D vector [FX, FY, FZ], got shape {self.force.shape}")


@dataclass
class NodalMoment:
    """
    Moment applied directly at a node (in global coordinates).
    
    Attributes:
        node_id: Target node ID
        moment: Moment vector [MX, MY, MZ] in global coordinates
    """
    node_id: str
    moment: np.ndarray
    
    def __post_init__(self) -> None:
        """Validate the nodal moment."""
        if not isinstance(self.moment, np.ndarray):
            self.moment = np.array(self.moment, dtype=float)
        if self.moment.shape != (3,):
            raise ValueError(f"NodalMoment must be a 3D vector [MX, MY, MZ], got shape {self.moment.shape}")


@dataclass
class MemberPointForce:
    """
    Point force applied along a member (in local or global coordinates).
    
    Attributes:
        member_id: Target member ID
        position: Distance from start node (absolute) or fraction (relative, 0 to 1)
        force: Force vector [Fx, Fy, Fz]
        coords: Coordinate system - "local" or "global"
        position_type: "absolute" (distance from start) or "relative" (fraction 0-1)
    """
    member_id: str
    position: float
    force: np.ndarray
    coords: str = "local"
    position_type: str = "absolute"
    
    def __post_init__(self) -> None:
        """Validate the member point force."""
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)
        if self.force.shape != (3,):
            raise ValueError(f"MemberPointForce must be a 3D vector, got shape {self.force.shape}")
        
        if self.coords not in ("local", "global"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        
        if self.position_type not in ("absolute", "relative"):
            raise ValueError(f"position_type must be 'absolute' or 'relative', got '{self.position_type}'")
        
        if self.position_type == "relative" and not (0 <= self.position <= 1):
            raise ValueError(f"Relative position must be in range [0, 1], got {self.position}")
        
        if self.position < 0:
            raise ValueError(f"Position must be non-negative, got {self.position}")


@dataclass
class MemberDistributedForce:
    """
    Distributed force along a member segment.
    
    Attributes:
        member_id: Target member ID
        start_position: Distance from start node (start of distributed load)
        end_position: Distance from start node (end of distributed load)
        start_force: Force per unit length [wx, wy, wz] at start position
        end_force: Force per unit length [wx, wy, wz] at end position
        coords: Coordinate system - "local" or "global"
    """
    member_id: str
    start_position: float
    end_position: float
    start_force: np.ndarray
    end_force: np.ndarray
    coords: str = "local"
    
    def __post_init__(self) -> None:
        """Validate the member distributed force."""
        if not isinstance(self.start_force, np.ndarray):
            self.start_force = np.array(self.start_force, dtype=float)
        if self.start_force.shape != (3,):
            raise ValueError(f"start_force must be a 3D vector, got shape {self.start_force.shape}")
        
        if not isinstance(self.end_force, np.ndarray):
            self.end_force = np.array(self.end_force, dtype=float)
        if self.end_force.shape != (3,):
            raise ValueError(f"end_force must be a 3D vector, got shape {self.end_force.shape}")
        
        if self.coords not in ("local", "global"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        
        if self.start_position < 0:
            raise ValueError(f"start_position must be non-negative, got {self.start_position}")
        
        # Allow end_position = -1 as a special marker for "full length"
        # LoadedFrame will replace -1 with actual member length during analysis
        if self.end_position != -1.0 and self.end_position < self.start_position:
            raise ValueError(
                f"end_position ({self.end_position}) must be >= start_position ({self.start_position}) or -1 (full length)"
            )


@dataclass
class FrameLoadCase:
    """
    Loads applied to a frame structure.
    
    Attributes:
        name: Load case identifier
        nodal_forces: List of forces applied at nodes
        nodal_moments: List of moments applied at nodes
        member_point_forces: List of point forces applied to members
        member_distributed_forces: List of distributed forces applied to members
    """
    name: str
    nodal_forces: List[NodalForce] = field(default_factory=list)
    nodal_moments: List[NodalMoment] = field(default_factory=list)
    member_point_forces: List[MemberPointForce] = field(default_factory=list)
    member_distributed_forces: List[MemberDistributedForce] = field(default_factory=list)
    
    def add_nodal_force(self, node_id: str, force: np.ndarray) -> None:
        """
        Add a force at a node.
        
        Args:
            node_id: Target node ID
            force: Force vector [FX, FY, FZ] in global coordinates
        """
        self.nodal_forces.append(NodalForce(node_id=node_id, force=force))
    
    def add_nodal_moment(self, node_id: str, moment: np.ndarray) -> None:
        """
        Add a moment at a node.
        
        Args:
            node_id: Target node ID
            moment: Moment vector [MX, MY, MZ] in global coordinates
        """
        self.nodal_moments.append(NodalMoment(node_id=node_id, moment=moment))
    
    def add_member_point_force(
        self,
        member_id: str,
        position: float,
        force: np.ndarray,
        coords: str = "local",
        position_type: str = "absolute"
    ) -> None:
        """
        Add a point force along a member.
        
        Args:
            member_id: Target member ID
            position: Distance from start node (if absolute) or fraction 0-1 (if relative)
            force: Force vector [Fx, Fy, Fz]
            coords: "local" (member coordinates) or "global"
            position_type: "absolute" (distance) or "relative" (fraction)
        """
        self.member_point_forces.append(
            MemberPointForce(
                member_id=member_id,
                position=position,
                force=force,
                coords=coords,
                position_type=position_type
            )
        )
    
    def add_member_distributed_force(
        self,
        member_id: str,
        start_position: float,
        end_position: float,
        start_force: np.ndarray,
        end_force: np.ndarray,
        coords: str = "local"
    ) -> None:
        """
        Add a distributed force along a member segment.
        
        Args:
            member_id: Target member ID
            start_position: Distance from start node (start of load)
            end_position: Distance from start node (end of load)
            start_force: Force per unit length [wx, wy, wz] at start
            end_force: Force per unit length [wx, wy, wz] at end
            coords: "local" (member coordinates) or "global"
        """
        self.member_distributed_forces.append(
            MemberDistributedForce(
                member_id=member_id,
                start_position=start_position,
                end_position=end_position,
                start_force=start_force,
                end_force=end_force,
                coords=coords
            )
        )
    
    def add_member_uniform_force(
        self,
        member_id: str,
        force: np.ndarray,
        coords: str = "local"
    ) -> None:
        """
        Add a uniform distributed force over the entire member length.
        
        Note: The actual member length is not known here, so we use placeholders.
        The LoadedFrame will adjust start_position=0 and end_position=L during analysis.
        
        Args:
            member_id: Target member ID
            force: Force per unit length [wx, wy, wz]
            coords: "local" (member coordinates) or "global"
        """
        # Use a special marker: start=0, end=-1 means "full length"
        # LoadedFrame will replace -1 with actual member length
        self.member_distributed_forces.append(
            MemberDistributedForce(
                member_id=member_id,
                start_position=0.0,
                end_position=-1.0,  # Marker for "full length"
                start_force=force,
                end_force=force,
                coords=coords
            )
        )
    
    def __repr__(self) -> str:
        """String representation of the load case."""
        n_nodal = len(self.nodal_forces) + len(self.nodal_moments)
        n_member = len(self.member_point_forces) + len(self.member_distributed_forces)
        return f"FrameLoadCase('{self.name}', {n_nodal} nodal loads, {n_member} member loads)"
