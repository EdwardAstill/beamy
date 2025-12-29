from __future__ import annotations

from typing import Tuple

import numpy as np

from beamy.analysis.transformations import XYZ, frame_transformation, length_and_direction, transform_load, transform_stiffness


def material_stiffness_local(
    e_modulus: float,
    g_modulus: float,
    area: float,
    iy: float,
    iz: float,
    j: float,
    length: float,
) -> np.ndarray:
    """
    3D Euler-Bernoulli space-frame local stiffness (12x12).

    Local DOF per node: [u, v, w, rx, ry, rz]
    """
    if length <= 0.0:
        raise ValueError("length must be positive")
    if area <= 0.0 or e_modulus <= 0.0:
        raise ValueError("area and e_modulus must be positive")
    if g_modulus <= 0.0:
        raise ValueError("g_modulus must be positive")
    if iy <= 0.0 or iz <= 0.0 or j <= 0.0:
        raise ValueError("iy, iz, j must be positive for beam elements")

    l = length
    k = np.zeros((12, 12), dtype=float)

    # axial
    ka = e_modulus * area / l
    k[0, 0] = ka
    k[0, 6] = -ka
    k[6, 0] = -ka
    k[6, 6] = ka

    # torsion
    kt = g_modulus * j / l
    k[3, 3] = kt
    k[3, 9] = -kt
    k[9, 3] = -kt
    k[9, 9] = kt

    # bending about local z (v, rz) uses Iz
    ez = e_modulus * iz
    k1 = 12.0 * ez / (l**3)
    k2 = 6.0 * ez / (l**2)
    k3 = 4.0 * ez / l
    k4 = 2.0 * ez / l
    # indices: v_i=1, rz_i=5, v_j=7, rz_j=11
    _add_4x4(k, (1, 5, 7, 11), np.array([[k1, k2, -k1, k2], [k2, k3, -k2, k4], [-k1, -k2, k1, -k2], [k2, k4, -k2, k3]]))

    # bending about local y (w, ry) uses Iy (note sign pattern)
    ey = e_modulus * iy
    k1 = 12.0 * ey / (l**3)
    k2 = 6.0 * ey / (l**2)
    k3 = 4.0 * ey / l
    k4 = 2.0 * ey / l
    # indices: w_i=2, ry_i=4, w_j=8, ry_j=10
    _add_4x4(
        k,
        (2, 4, 8, 10),
        np.array([[k1, -k2, -k1, -k2], [-k2, k3, k2, k4], [-k1, k2, k1, k2], [-k2, k4, k2, k3]]),
    )

    return k


def geometric_stiffness_local(axial_force_value: float, length: float) -> np.ndarray:
    """
    Simplified geometric stiffness (initial-stress) for a 3D frame element.

    Adds geometric stiffness in the bending DOFs only using the classic
    2D beam-column initial-stress matrix applied to both bending planes.

    Uses N>0 tension, N<0 compression.
    """
    if length <= 0.0:
        raise ValueError("length must be positive")

    n = float(axial_force_value)
    l = float(length)
    kg = np.zeros((12, 12), dtype=float)
    if abs(n) <= 0.0:
        return kg

    factor = n / (30.0 * l)
    base = np.array(
        [
            [36.0, 3.0 * l, -36.0, 3.0 * l],
            [3.0 * l, 4.0 * l * l, -3.0 * l, -1.0 * l * l],
            [-36.0, -3.0 * l, 36.0, -3.0 * l],
            [3.0 * l, -1.0 * l * l, -3.0 * l, 4.0 * l * l],
        ],
        dtype=float,
    )

    # plane v-rz
    _add_4x4(kg, (1, 5, 7, 11), factor * base)
    # plane w-ry
    _add_4x4(kg, (2, 4, 8, 10), factor * base)

    return kg


def equivalent_nodal_loads_local_uniform(w_local: Tuple[float, float, float], length: float) -> np.ndarray:
    """
    Consistent nodal load vector (local) for uniform distributed load along the member.

    w_local = (wx, wy, wz) force/length in local axes.
    """
    wx, wy, wz = float(w_local[0]), float(w_local[1]), float(w_local[2])
    l = float(length)
    f = np.zeros(12, dtype=float)

    # axial
    f[0] += wx * l / 2.0
    f[6] += wx * l / 2.0

    # transverse local y (couples into Mz)
    f[1] += wy * l / 2.0
    f[7] += wy * l / 2.0
    f[5] += wy * (l**2) / 12.0
    f[11] -= wy * (l**2) / 12.0

    # transverse local z (couples into My)
    f[2] += wz * l / 2.0
    f[8] += wz * l / 2.0
    f[4] -= wz * (l**2) / 12.0
    f[10] += wz * (l**2) / 12.0

    return f


def material_stiffness_global(p0: XYZ, p1: XYZ, e_modulus: float, g_modulus: float, area: float, iy: float, iz: float, j: float) -> np.ndarray:
    length, _ = length_and_direction(p0, p1)
    k_local = material_stiffness_local(e_modulus, g_modulus, area, iy, iz, j, length)
    t = frame_transformation(p0, p1)
    return transform_stiffness(k_local, t)


def geometric_stiffness_global(p0: XYZ, p1: XYZ, axial_force_value: float) -> np.ndarray:
    length, _ = length_and_direction(p0, p1)
    kg_local = geometric_stiffness_local(axial_force_value, length)
    t = frame_transformation(p0, p1)
    return transform_stiffness(kg_local, t)


def equivalent_nodal_loads_global_uniform(p0: XYZ, p1: XYZ, w_local: Tuple[float, float, float]) -> np.ndarray:
    length, _ = length_and_direction(p0, p1)
    f_local = equivalent_nodal_loads_local_uniform(w_local, length)
    t = frame_transformation(p0, p1)
    return transform_load(f_local, t)


def _add_4x4(k: np.ndarray, dof: Tuple[int, int, int, int], sub: np.ndarray) -> None:
    a, b, c, d = dof
    idx = (a, b, c, d)
    for i in range(4):
        for j in range(4):
            k[idx[i], idx[j]] += float(sub[i, j])

