from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union

from beamy.model.member import MemberKind

XYZ = Tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class AnalysisNode:
    id: int
    source_id: Union[int, str]
    xyz: XYZ


@dataclass(frozen=True, slots=True)
class Element:
    id: int
    node_i: int
    node_j: int
    parent_member_id: Union[int, str]
    kind: MemberKind
    area: float
    youngs_modulus: float
    shear_modulus: float
    iy: float
    iz: float
    j: float


@dataclass
class AnalysisModel:
    nodes: List[AnalysisNode] = field(default_factory=list)
    elements: List[Element] = field(default_factory=list)
    member_to_elements: Dict[Union[int, str], List[int]] = field(default_factory=dict)
    node_lookup: Dict[str, int] = field(default_factory=dict)

    def add_node(self, source_id: Union[int, str], xyz: XYZ) -> int:
        node_id = len(self.nodes)
        self.nodes.append(AnalysisNode(node_id, source_id, xyz))
        self.node_lookup[str(source_id)] = node_id
        return node_id

    def add_element(
        self,
        node_i: int,
        node_j: int,
        parent_member_id: Union[int, str],
        kind: MemberKind,
        area: float,
        youngs_modulus: float,
        shear_modulus: float,
        iy: float,
        iz: float,
        j: float,
    ) -> int:
        element_id = len(self.elements)
        self.elements.append(
            Element(
                id=element_id,
                node_i=node_i,
                node_j=node_j,
                parent_member_id=parent_member_id,
                kind=kind,
                area=area,
                youngs_modulus=youngs_modulus,
                shear_modulus=shear_modulus,
                iy=iy,
                iz=iz,
                j=j,
            )
        )
        if parent_member_id not in self.member_to_elements:
            self.member_to_elements[parent_member_id] = []
        self.member_to_elements[parent_member_id].append(element_id)
        return element_id

