from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional
import numpy as np
from .node import Node
from .member import Member
from .mpc import MPC

@dataclass
class Frame:
    """A collection of nodes and members forming a 3D structural system."""
    members: List[Member]
    mpcs: List[MPC] = field(default_factory=list)
    nodes: Dict[str, Node] = field(init=False, default_factory=dict)
    
    def __post_init__(self) -> None:
        raise NotImplementedError("Use Frame.from_nodes_and_members(nodes, members)")
    
    @classmethod
    def from_nodes_and_members(cls, nodes: List[Node], members: List[Member], mpcs: Optional[List[MPC]] = None) -> Frame:
        frame = cls.__new__(cls)
        frame.members = members
        frame.mpcs = mpcs or []
        frame.nodes = {}
        for n in nodes:
            if n.id in frame.nodes: raise ValueError(f"Duplicate node ID: {n.id}")
            frame.nodes[n.id] = n
        frame._validate_and_setup_members()
        frame._validate_supports()
        return frame
    
    def _validate_and_setup_members(self) -> None:
        mids = set()
        for m in self.members:
            if m.id in mids: raise ValueError(f"Duplicate member ID: {m.id}")
            mids.add(m.id)
            if m.start_node_id not in self.nodes: raise ValueError(f"Member {m.id}: start_node {m.start_node_id} not found")
            if m.end_node_id not in self.nodes: raise ValueError(f"Member {m.id}: end_node {m.end_node_id} not found")
            m.set_nodes(self.nodes[m.start_node_id], self.nodes[m.end_node_id])
            if m.length < 1e-10: raise ValueError(f"Member {m.id}: zero-length")
            if m.id not in self.nodes[m.start_node_id].connected_members: self.nodes[m.start_node_id].connected_members.append(m.id)
            if m.id not in self.nodes[m.end_node_id].connected_members: self.nodes[m.end_node_id].connected_members.append(m.id)
    
    def _validate_supports(self) -> None:
        constrained_dofs: Set[Tuple[str, int]] = set()

        for node in self.nodes.values():
            if node.support:
                for d in range(6):
                    if node.support[d] == "1":
                        constrained_dofs.add((node.id, d))

        for m in self.members:
            if m.constraints:
                for d in range(6):
                    if m.constraints[d] == "1":
                        constrained_dofs.add((m.start_node_id, d))
                    if m.constraints[6 + d] == "1":
                        constrained_dofs.add((m.end_node_id, d))

        if not constrained_dofs:
            raise ValueError("No supports found (neither Node.support nor Member.constraints)")
        if len(constrained_dofs) < 6:
            raise ValueError("Insufficient supports for 3D stability (< 6 constrained DOFs)")
    
    @property
    def node_positions(self) -> Dict[str, np.ndarray]: return {nid: n.position for nid, n in self.nodes.items()}
    @property
    def supported_nodes(self) -> List[Node]: return [n for n in self.nodes.values() if n.support]
    @property
    def member_lengths(self) -> Dict[str, float]: return {m.id: m.length for m in self.members}
    
    def get_member(self, mid: str) -> Member:
        for m in self.members:
            if m.id == mid: return m
        raise KeyError(f"Member '{mid}' not found")
    
    def get_node(self, nid: str) -> Node:
        if nid not in self.nodes: raise KeyError(f"Node '{nid}' not found")
        return self.nodes[nid]
    
    def members_at_node(self, nid: str) -> List[Member]:
        return [self.get_member(mid) for mid in self.get_node(nid).connected_members]
    
    def __repr__(self) -> str: return f"Frame({len(self.nodes)} nodes, {len(self.members)} members)"
