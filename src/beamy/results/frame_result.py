from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple, Union

from beamy.results.member_result import MemberResult

Disp6 = Tuple[float, float, float, float, float, float]


@dataclass
class FrameResult:
    loadcase: str
    node_displacements: Dict[str, Disp6] = field(default_factory=dict)
    reactions: Dict[str, Disp6] = field(default_factory=dict)
    member_results: Dict[Union[int, str], MemberResult] = field(default_factory=dict)

    def displacement_at(self, node_id: Union[int, str]) -> Disp6:
        key = str(node_id)
        if key not in self.node_displacements:
            raise KeyError(f"displacement for node {node_id} not found")
        return self.node_displacements[key]

    def reaction_at(self, node_id: Union[int, str]) -> Disp6:
        key = str(node_id)
        if key not in self.reactions:
            raise KeyError(f"reaction for node {node_id} not found")
        return self.reactions[key]

    def member(self, member_id: Union[int, str]) -> MemberResult:
        if member_id not in self.member_results:
            raise KeyError(f"member result for {member_id} not found")
        return self.member_results[member_id]

