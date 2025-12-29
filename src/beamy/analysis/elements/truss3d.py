from __future__ import annotations

from typing import Iterable, Tuple

import numpy as np

from beamy.analysis.transformations import length_and_direction, truss_transformation

XYZ = Tuple[float, float, float]


def local_stiffness_6x6(e_modulus: float, area: float, length: float) -> np.ndarray:
    if length <= 0.0:
        raise ValueError("element length must be positive")
    k = e_modulus * area / length
    mat = np.zeros((6, 6), dtype=float)
    mat[0, 0] = k
    mat[0, 3] = -k
    mat[3, 0] = -k
    mat[3, 3] = k
    return mat


def global_stiffness_6x6(p0: XYZ, p1: XYZ, e_modulus: float, area: float) -> np.ndarray:
    length, _ = length_and_direction(p0, p1)
    k_local = local_stiffness_6x6(e_modulus, area, length)
    t = truss_transformation(p0, p1)
    return t.T @ k_local @ t


def global_stiffness_12x12(p0: XYZ, p1: XYZ, e_modulus: float, area: float) -> np.ndarray:
    """
    Embed translational truss stiffness into a 12x12 (6 dof/node) matrix.
    DOF order per node: [UX, UY, UZ, RX, RY, RZ].
    """
    k6 = global_stiffness_6x6(p0, p1, e_modulus, area)
    k12 = np.zeros((12, 12), dtype=float)
    map_idx = [0, 1, 2, 6, 7, 8]
    for a in range(6):
        for b in range(6):
            k12[map_idx[a], map_idx[b]] = k6[a, b]
    return k12


def geometric_stiffness_12x12(p0: XYZ, p1: XYZ, axial_force_value: float) -> np.ndarray:
    """
    Initial-stress geometric stiffness for truss (embedded 12x12).

    Uses N>0 tension, N<0 compression.
    """
    length, direction = length_and_direction(p0, p1)
    l = np.array(direction, dtype=float).reshape(3, 1)
    ident = np.eye(3, dtype=float)
    p = ident - (l @ l.T)
    kg_block = (axial_force_value / length) * p
    kg6 = np.zeros((6, 6), dtype=float)
    kg6[0:3, 0:3] = kg_block
    kg6[0:3, 3:6] = -kg_block
    kg6[3:6, 0:3] = -kg_block
    kg6[3:6, 3:6] = kg_block

    kg12 = np.zeros((12, 12), dtype=float)
    map_idx = [0, 1, 2, 6, 7, 8]
    for a in range(6):
        for b in range(6):
            kg12[map_idx[a], map_idx[b]] = kg6[a, b]
    return kg12


def axial_force(local_displacements: Iterable[float], e_modulus: float, area: float, length: float) -> float:
    disp = tuple(float(value) for value in local_displacements)
    if len(disp) != 6:
        raise ValueError("local_displacements must have 6 components")
    du = disp[3] - disp[0]
    return (e_modulus * area / length) * du

