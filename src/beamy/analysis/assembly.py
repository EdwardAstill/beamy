from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from beamy.analysis.elements import element_matrices
from beamy.analysis.model import AnalysisModel
from beamy.analysis.transformations import local_rotation_matrix
from beamy.analysis.elements.beam3d import equivalent_nodal_loads_global_uniform
from beamy.loads.loadcase import LoadCase
from beamy.loads.member_loads import MemberDistributedLoad
from beamy.model.frame import Frame


@dataclass
class AssemblyResult:
    k_full: np.ndarray
    f_full: np.ndarray
    free_dofs: List[int]
    restrained_dofs: List[int]
    mask: List[bool]
    prescribed: Dict[int, float]


def assemble_global_system(
    model: AnalysisModel,
    frame: Frame,
    loadcase: LoadCase,
    element_active: Optional[Sequence[bool]] = None,
    geometric_axial_forces: Optional[Sequence[float]] = None,
    include_geometric: bool = False,
) -> AssemblyResult:
    total_dofs = len(model.nodes) * 6
    mask = [False for _ in range(total_dofs)]
    prescribed: Dict[int, float] = {}

    # If there are no beam elements, rotations have no stiffness; constrain them to avoid singular systems.
    has_beam = any(element.kind == "beam" for element in model.elements)
    if not has_beam:
        for node_id in range(len(model.nodes)):
            base = node_id * 6
            mask[base + 3] = True
            mask[base + 4] = True
            mask[base + 5] = True

    # Apply supports
    for support in loadcase.supports:
        node_key = str(support.target)
        if node_key not in model.node_lookup:
            raise ValueError(f"support node {support.target} missing in analysis model")
        node_idx = model.node_lookup[node_key]
        base = node_idx * 6
        for local_dof in range(6):
            if support.restrained_dofs[local_dof]:
                mask[base + local_dof] = True

    # Prescribed displacements
    for disp in loadcase.prescribed_displacements:
        node_key = str(disp.target)
        if node_key not in model.node_lookup:
            raise ValueError(f"prescribed displacement node {disp.target} missing in analysis model")
        node_idx = model.node_lookup[node_key]
        base = node_idx * 6
        for local_dof in range(6):
            dof_idx = base + local_dof
            mask[dof_idx] = True
            prescribed[dof_idx] = float(disp.values[local_dof])

    k_full = np.zeros((total_dofs, total_dofs), dtype=float)
    f_full = np.zeros(total_dofs, dtype=float)

    # Assemble element stiffness
    active = list(element_active) if element_active is not None else [True for _ in model.elements]
    if len(active) != len(model.elements):
        raise ValueError("element_active length must match number of elements")

    axial_forces = list(geometric_axial_forces) if geometric_axial_forces is not None else [0.0 for _ in model.elements]
    if geometric_axial_forces is not None and len(axial_forces) != len(model.elements):
        raise ValueError("geometric_axial_forces length must match number of elements")

    for element in model.elements:
        if not active[element.id]:
            continue
        node_i = model.nodes[element.node_i]
        node_j = model.nodes[element.node_j]
        mats = element_matrices(
            element,
            node_i.xyz,
            node_j.xyz,
            axial_force_value=axial_forces[element.id],
            include_geometric=include_geometric,
        )

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
        for a in range(12):
            for b in range(12):
                k_full[dof_indices[a], dof_indices[b]] += float(mats.k_global[a, b])
        for a in range(12):
            f_full[dof_indices[a]] += float(mats.f_equiv_global[a])

    # Assemble nodal loads
    for load in loadcase.nodal_loads:
        node_key = str(load.node_id)
        if node_key not in model.node_lookup:
            raise ValueError(f"nodal load references missing node {load.node_id}")
        idx = model.node_lookup[node_key] * 6
        f_full[idx + 0] += float(load.forces[0])
        f_full[idx + 1] += float(load.forces[1])
        f_full[idx + 2] += float(load.forces[2])
        f_full[idx + 3] += float(load.moments[0])
        f_full[idx + 4] += float(load.moments[1])
        f_full[idx + 5] += float(load.moments[2])

    # Assemble member loads (MVP: full-length uniform distributed load only)
    for mload in loadcase.member_loads:
        if isinstance(mload, MemberDistributedLoad):
            if mload.s0 != 0.0 or mload.s1 != 1.0:
                raise NotImplementedError("member distributed loads supported only for s0=0, s1=1 in MVP")
            if mload.w0 != mload.w1:
                raise NotImplementedError("member distributed loads supported only for uniform w0==w1 in MVP")
            if mload.member_id not in frame.members:
                raise ValueError(f"member load references missing member {mload.member_id}")
            member = frame.members[mload.member_id]
            if mload.member_id not in model.member_to_elements:
                raise ValueError(f"member {mload.member_id} missing from analysis model mapping")
            element_ids = model.member_to_elements[mload.member_id]
            if len(element_ids) != 1:
                raise NotImplementedError("member load supported only when member maps to exactly one element in MVP")
            element = model.elements[element_ids[0]]
            if element.kind != "beam":
                raise NotImplementedError("member distributed loads implemented only for beam elements in MVP")

            node_i = model.nodes[element.node_i]
            node_j = model.nodes[element.node_j]

            w = (float(mload.w0[0]), float(mload.w0[1]), float(mload.w0[2]))
            if mload.coord_sys == "global":
                rot = local_rotation_matrix(node_i.xyz, node_j.xyz)
                w_local_arr = rot @ np.array(w, dtype=float)
                w_local = (float(w_local_arr[0]), float(w_local_arr[1]), float(w_local_arr[2]))
            else:
                w_local = w

            f_equiv = equivalent_nodal_loads_global_uniform(node_i.xyz, node_j.xyz, w_local)
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
            for a in range(12):
                f_full[dof_indices[a]] += float(f_equiv[a])
        else:
            raise NotImplementedError("only MemberDistributedLoad is supported in MVP member load assembly")

    free_dofs = [i for i, restrained in enumerate(mask) if not restrained]
    restrained_dofs = [i for i, restrained in enumerate(mask) if restrained]

    return AssemblyResult(
        k_full=k_full,
        f_full=f_full,
        free_dofs=free_dofs,
        restrained_dofs=restrained_dofs,
        mask=mask,
        prescribed=prescribed,
    )

