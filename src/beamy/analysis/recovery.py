from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from beamy.analysis.assembly import AssemblyResult
from beamy.analysis.model import AnalysisModel
from beamy.analysis.elements.beam3d import material_stiffness_local
from beamy.analysis.transformations import frame_transformation, length_and_direction
from beamy.results.frame_result import FrameResult
from beamy.results.member_result import MemberResult

Disp6 = Tuple[float, float, float, float, float, float]


def recover_results(model: AnalysisModel, assembly: AssemblyResult, displacements: np.ndarray, loadcase_name: str) -> FrameResult:
    node_displacements = _collect_node_displacements(model, displacements)
    reactions = _collect_reactions(model, assembly, displacements)
    member_results = _recover_members(model, displacements)

    return FrameResult(
        loadcase=loadcase_name,
        node_displacements=node_displacements,
        reactions=reactions,
        member_results=member_results,
    )


def _collect_node_displacements(model: AnalysisModel, displacements: np.ndarray) -> Dict[str, Disp6]:
    result: Dict[str, Disp6] = {}
    for node in model.nodes:
        base = node.id * 6
        disp = (
            float(displacements[base + 0]),
            float(displacements[base + 1]),
            float(displacements[base + 2]),
            float(displacements[base + 3]),
            float(displacements[base + 4]),
            float(displacements[base + 5]),
        )
        result[str(node.source_id)] = disp
    return result


def _collect_reactions(model: AnalysisModel, assembly: AssemblyResult, displacements: np.ndarray) -> Dict[str, Disp6]:
    full_reactions = assembly.k_full @ displacements - assembly.f_full
    reactions: Dict[str, Disp6] = {}
    for node in model.nodes:
        base = node.id * 6
        r = (
            float(full_reactions[base + 0]),
            float(full_reactions[base + 1]),
            float(full_reactions[base + 2]),
            float(full_reactions[base + 3]),
            float(full_reactions[base + 4]),
            float(full_reactions[base + 5]),
        )
        if any(abs(value) > 1e-9 for value in r):
            reactions[str(node.source_id)] = r
    return reactions


def _recover_members(model: AnalysisModel, displacements: np.ndarray) -> Dict[Union[int, str], MemberResult]:
    results: Dict[Union[int, str], MemberResult] = {}
    axial_forces = element_axial_forces(model, displacements)
    for element in model.elements:
        end_forces_local = _element_end_forces_local(element, model, displacements)
        results[element.parent_member_id] = MemberResult(
            member_id=element.parent_member_id,
            kind=element.kind,
            axial_force=float(axial_forces[element.id]),
            end_forces_local=end_forces_local,
        )
    return results


def element_axial_forces(model: AnalysisModel, displacements: np.ndarray) -> List[float]:
    """
    Compute axial force N for each element using local axial displacements.

    Convention:
    - N > 0 tension
    - N < 0 compression
    """
    forces: List[float] = []
    for element in model.elements:
        node_i = model.nodes[element.node_i]
        node_j = model.nodes[element.node_j]
        length, _ = length_and_direction(node_i.xyz, node_j.xyz)

        dof_indices = [
            element.node_i * 6 + 0,
            element.node_i * 6 + 1,
            element.node_i * 6 + 2,
            element.node_i * 6 + 3,
            element.node_i * 6 + 4,
            element.node_i * 6 + 5,
            element.node_j * 6 + 0,
            element.node_j * 6 + 1,
            element.node_j * 6 + 2,
            element.node_j * 6 + 3,
            element.node_j * 6 + 4,
            element.node_j * 6 + 5,
        ]
        u_global = np.array([float(displacements[i]) for i in dof_indices], dtype=float)
        t = frame_transformation(node_i.xyz, node_j.xyz)
        u_local = t @ u_global

        du = float(u_local[6] - u_local[0])
        n = (element.youngs_modulus * element.area / length) * du
        forces.append(float(n))
    return forces


def _element_end_forces_local(element, model: AnalysisModel, displacements: np.ndarray) -> Optional[Tuple[float, ...]]:
    """
    Compute basic local end force vector for beams; for truss/cable return axial-only vector (6 dof embed).
    """
    node_i = model.nodes[element.node_i]
    node_j = model.nodes[element.node_j]

    dof_indices = [
        element.node_i * 6 + 0,
        element.node_i * 6 + 1,
        element.node_i * 6 + 2,
        element.node_i * 6 + 3,
        element.node_i * 6 + 4,
        element.node_i * 6 + 5,
        element.node_j * 6 + 0,
        element.node_j * 6 + 1,
        element.node_j * 6 + 2,
        element.node_j * 6 + 3,
        element.node_j * 6 + 4,
        element.node_j * 6 + 5,
    ]
    u_global = np.array([float(displacements[i]) for i in dof_indices], dtype=float)
    t = frame_transformation(node_i.xyz, node_j.xyz)
    u_local = t @ u_global
    length, _ = length_and_direction(node_i.xyz, node_j.xyz)

    if element.kind == "beam":
        k_local = material_stiffness_local(
            element.youngs_modulus,
            element.shear_modulus,
            element.area,
            element.iy,
            element.iz,
            element.j,
            length,
        )
        f_local = k_local @ u_local
        return tuple(float(x) for x in f_local.tolist())

    # truss/cable: axial only in local x translation (u)
    n = (element.youngs_modulus * element.area / length) * float(u_local[6] - u_local[0])
    # return a sparse 12-vector consistent with frame local dof ordering
    f = [0.0 for _ in range(12)]
    f[0] = -float(n)
    f[6] = float(n)
    return tuple(f)

