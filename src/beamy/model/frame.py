from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Union

from beamy.model.connections import Connection
from beamy.model.member import Member
from beamy.model.node import Node, NodeId


@dataclass
class Frame:
    """
    User-facing structural model.

    Holds nodes, members, and connection references. No solver state is stored here.
    """

    nodes: Dict[NodeId, Node] = field(default_factory=dict)
    members: Dict[Union[int, str], Member] = field(default_factory=dict)
    connections: List[Connection] = field(default_factory=list)
    merge_tol: float = 1e-6

    def add_node(self, node: Node) -> None:
        if node.id in self.nodes:
            raise ValueError(f"node {node.id} already exists")
        self.nodes[node.id] = node

    def add_member(self, member: Member) -> None:
        if member.id in self.members:
            raise ValueError(f"member {member.id} already exists")
        if member.end_i not in self.nodes or member.end_j not in self.nodes:
            raise ValueError("member endpoints must reference existing nodes")
        self.members[member.id] = member

    def remove_node(self, node_id: NodeId) -> None:
        if node_id not in self.nodes:
            raise KeyError(f"node {node_id} not found")
        if any(member.end_i == node_id or member.end_j == node_id for member in self.members.values()):
            raise ValueError("cannot remove node while members reference it")
        self.nodes.pop(node_id)

    def remove_member(self, member_id: Union[int, str]) -> None:
        if member_id not in self.members:
            raise KeyError(f"member {member_id} not found")
        self.members.pop(member_id)

    def validate(self) -> None:
        self._validate_nodes()
        self._validate_members()

    def analyze(self, loadcase: "LoadCase", settings: Optional["AnalysisSettings"] = None) -> "FrameResult":
        from beamy.analysis import run_analysis
        from beamy.analysis.settings import AnalysisSettings

        chosen = settings or AnalysisSettings()
        return run_analysis(self, loadcase, chosen)

    def plot(self, save_path: Optional[str] = None) -> None:
        from beamy.viz.plot_frame import plot_frame_undeformed

        plot_frame_undeformed(self, save_path=save_path)

    def _validate_nodes(self) -> None:
        if not self.nodes:
            raise ValueError("frame must contain at least one node")
        for node in self.nodes.values():
            if not isinstance(node, Node):
                raise TypeError("all nodes must be Node instances")

    def _validate_members(self) -> None:
        if not self.members:
            raise ValueError("frame must contain at least one member")
        for member in self.members.values():
            if not isinstance(member, Member):
                raise TypeError("all members must be Member instances")
            if member.end_i not in self.nodes or member.end_j not in self.nodes:
                raise ValueError(f"member {member.id} references missing nodes")

