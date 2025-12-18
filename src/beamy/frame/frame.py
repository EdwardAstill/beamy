# frame.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict
import numpy as np

from .node import Node
from .member import Member


@dataclass
class Frame:
    """
    A collection of nodes and members forming a 3D structural system.
    
    Attributes:
        members: List of Member objects defining the frame geometry
        nodes: Dictionary mapping node IDs to Node objects (auto-built from members)
    """
    members: List[Member]
    nodes: Dict[str, Node] = field(init=False, default_factory=dict)
    
    def __post_init__(self) -> None:
        """Build node dictionary and validate frame after initialization."""
        # This constructor expects nodes to be built from member references
        # For explicit node+member construction, use from_nodes_and_members()
        raise NotImplementedError(
            "Frame() constructor requires explicit nodes. "
            "Use Frame.from_nodes_and_members(nodes, members) instead."
        )
    
    @classmethod
    def from_nodes_and_members(cls, nodes: List[Node], members: List[Member]) -> Frame:
        """
        Create a Frame from explicit node and member lists.
        
        Args:
            nodes: List of Node objects
            members: List of Member objects
            
        Returns:
            Frame object with validated geometry
            
        Raises:
            ValueError: If validation fails (duplicate IDs, missing references, etc.)
        """
        # Create frame instance (bypass __post_init__ validation)
        frame = cls.__new__(cls)
        frame.members = members
        
        # Build nodes dictionary
        frame.nodes = {}
        for node in nodes:
            if node.id in frame.nodes:
                raise ValueError(f"Duplicate node ID: {node.id}")
            frame.nodes[node.id] = node
        
        # Validate and setup members
        frame._validate_and_setup_members()
        
        # Validate supports
        frame._validate_supports()
        
        return frame
    
    def _validate_and_setup_members(self) -> None:
        """
        Validate member references and setup connections.
        
        Raises:
            ValueError: If member validation fails
        """
        member_ids = set()
        
        for member in self.members:
            # Check for duplicate member IDs
            if member.id in member_ids:
                raise ValueError(f"Duplicate member ID: {member.id}")
            member_ids.add(member.id)
            
            # Check that referenced nodes exist
            if member.start_node_id not in self.nodes:
                raise ValueError(
                    f"Member {member.id}: start_node_id '{member.start_node_id}' not found in nodes"
                )
            if member.end_node_id not in self.nodes:
                raise ValueError(
                    f"Member {member.id}: end_node_id '{member.end_node_id}' not found in nodes"
                )
            
            # Set node references in member
            start_node = self.nodes[member.start_node_id]
            end_node = self.nodes[member.end_node_id]
            member.set_nodes(start_node, end_node)
            
            # Check for zero-length members
            length = member.length
            if length < 1e-10:
                raise ValueError(
                    f"Member {member.id}: zero-length member "
                    f"(nodes {member.start_node_id} and {member.end_node_id} at same position)"
                )
            
            # Add member to node's connected_members list
            if member.id not in start_node.connected_members:
                start_node.connected_members.append(member.id)
            if member.id not in end_node.connected_members:
                end_node.connected_members.append(member.id)
    
    def _validate_supports(self) -> None:
        """
        Validate that frame has sufficient supports for stability.
        
        Basic check: at least one support exists and provides some constraint.
        More sophisticated stability checks would require eigenvalue analysis.
        
        Raises:
            ValueError: If no supports found
        """
        supported_nodes = [n for n in self.nodes.values() if n.support is not None]
        
        if not supported_nodes:
            raise ValueError("Frame has no supports. At least one node must have a support.")
        
        # Check that at least one DOF is constrained globally
        any_constrained = False
        for node in supported_nodes:
            if '1' in node.support:
                any_constrained = True
                break
        
        if not any_constrained:
            raise ValueError("Frame supports do not constrain any DOFs (all supports are '000000')")
        
        # Basic stability check: count total constrained DOFs
        # For a 3D structure, we need at least 6 constraints to prevent rigid body motion
        # (3 translations + 3 rotations)
        total_constraints = sum(
            node.support.count('1') 
            for node in supported_nodes
        )
        
        if total_constraints < 6:
            raise ValueError(
                f"Insufficient supports for stability. "
                f"Found {total_constraints} constrained DOFs, need at least 6 "
                f"(3 translations + 3 rotations) to prevent rigid body motion."
            )
    
    @property
    def node_positions(self) -> Dict[str, np.ndarray]:
        """Dictionary of node ID → position array."""
        return {node_id: node.position for node_id, node in self.nodes.items()}
    
    @property
    def supported_nodes(self) -> List[Node]:
        """List of nodes that have supports defined."""
        return [node for node in self.nodes.values() if node.support is not None]
    
    @property
    def member_lengths(self) -> Dict[str, float]:
        """Dictionary of member ID → length."""
        return {member.id: member.length for member in self.members}
    
    def get_member(self, member_id: str) -> Member:
        """
        Retrieve a member by ID.
        
        Args:
            member_id: Member identifier
            
        Returns:
            Member object
            
        Raises:
            KeyError: If member_id not found
        """
        for member in self.members:
            if member.id == member_id:
                return member
        raise KeyError(f"Member '{member_id}' not found in frame")
    
    def get_node(self, node_id: str) -> Node:
        """
        Retrieve a node by ID.
        
        Args:
            node_id: Node identifier
            
        Returns:
            Node object
            
        Raises:
            KeyError: If node_id not found
        """
        if node_id not in self.nodes:
            raise KeyError(f"Node '{node_id}' not found in frame")
        return self.nodes[node_id]
    
    def members_at_node(self, node_id: str) -> List[Member]:
        """
        Get all members connected to a node.
        
        Args:
            node_id: Node identifier
            
        Returns:
            List of Member objects connected to the node
        """
        node = self.get_node(node_id)
        return [self.get_member(mid) for mid in node.connected_members]
    
    def __repr__(self) -> str:
        """String representation of the frame."""
        return f"Frame({len(self.nodes)} nodes, {len(self.members)} members)"
