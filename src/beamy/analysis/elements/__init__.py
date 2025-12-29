from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from beamy.analysis.elements.beam3d import (
    geometric_stiffness_global as beam_geometric_stiffness,
    material_stiffness_global as beam_material_stiffness,
)
from beamy.analysis.elements.truss3d import geometric_stiffness_12x12 as truss_geometric_stiffness
from beamy.analysis.elements.truss3d import global_stiffness_12x12 as truss_material_stiffness
from beamy.analysis.model import Element


@dataclass(frozen=True)
class ElementMatrices:
    k_global: np.ndarray  # 12x12
    f_equiv_global: np.ndarray  # 12,


def element_matrices(
    element: Element,
    node_i_xyz,
    node_j_xyz,
    axial_force_value: Optional[float] = None,
    include_geometric: bool = False,
) -> ElementMatrices:
    """
    Return element stiffness (12x12) and equivalent nodal loads (12,) in global coords.
    """
    f_equiv = np.zeros(12, dtype=float)

    if element.kind in ("truss", "cable"):
        k = truss_material_stiffness(node_i_xyz, node_j_xyz, element.youngs_modulus, element.area)
        if include_geometric and axial_force_value is not None:
            k = k + truss_geometric_stiffness(node_i_xyz, node_j_xyz, axial_force_value)
        return ElementMatrices(k_global=k, f_equiv_global=f_equiv)

    if element.kind == "beam":
        k = beam_material_stiffness(
            node_i_xyz,
            node_j_xyz,
            element.youngs_modulus,
            element.shear_modulus,
            element.area,
            element.iy,
            element.iz,
            element.j,
        )
        if include_geometric and axial_force_value is not None:
            k = k + beam_geometric_stiffness(node_i_xyz, node_j_xyz, axial_force_value)
        return ElementMatrices(k_global=k, f_equiv_global=f_equiv)

    raise ValueError(f"unknown element kind {element.kind}")

