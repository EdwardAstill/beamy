from __future__ import annotations

from typing import Dict

from beamy.analysis.model import AnalysisModel
from beamy.model.frame import Frame
from beamy.loads.loadcase import LoadCase


def build_analysis_model(frame: Frame, loadcase: LoadCase) -> AnalysisModel:
    """
    Very simple meshing:
    - Creates one analysis node per frame node.
    - Creates one element per member (no intermediate stations yet).
    """

    frame.validate()
    loadcase.validate(frame)

    model = AnalysisModel()
    node_index: Dict[str, int] = {}

    for node_id, node in frame.nodes.items():
        node_index[str(node_id)] = model.add_node(node_id, node.xyz)

    for member in frame.members.values():
        i_idx = node_index[str(member.end_i)]
        j_idx = node_index[str(member.end_j)]
        model.add_element(
            node_i=i_idx,
            node_j=j_idx,
            parent_member_id=member.id,
            kind=member.kind,
            area=member.area,
            youngs_modulus=member.youngs_modulus,
            shear_modulus=member.shear_modulus,
            iy=member.iy,
            iz=member.iz,
            j=member.j,
        )

    return model

