from __future__ import annotations

from typing import Tuple

import numpy as np

XYZ = Tuple[float, float, float]
Matrix3 = np.ndarray
Matrix6 = np.ndarray
Matrix12 = np.ndarray


def unit_vector(p0: XYZ, p1: XYZ) -> XYZ:
    dx = float(p1[0] - p0[0])
    dy = float(p1[1] - p0[1])
    dz = float(p1[2] - p0[2])
    length = (dx * dx + dy * dy + dz * dz) ** 0.5
    if length <= 0.0:
        raise ValueError("zero length vector")
    return (dx / length, dy / length, dz / length)


def length_and_direction(p0: XYZ, p1: XYZ) -> Tuple[float, XYZ]:
    direction = unit_vector(p0, p1)
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    dz = p1[2] - p0[2]
    length = float((dx * dx + dy * dy + dz * dz) ** 0.5)
    return length, direction


def local_rotation_matrix(p0: XYZ, p1: XYZ) -> Matrix3:
    """
    Build a 3x3 rotation matrix where local x aligns with the member axis.

    The returned matrix R maps global vectors into local components:
    v_local = R @ v_global
    """

    x_local = unit_vector(p0, p1)
    up = (0.0, 0.0, 1.0)
    if abs(x_local[0]) < 1e-6 and abs(x_local[1]) < 1e-6:
        up = (0.0, 1.0, 0.0)

    y_local = _normalize(_cross(up, x_local))
    if _norm(y_local) < 1e-8:
        up = (0.0, 1.0, 0.0)
        y_local = _normalize(_cross(up, x_local))
    z_local = _cross(x_local, y_local)

    return np.array(
        [
            [x_local[0], x_local[1], x_local[2]],
            [y_local[0], y_local[1], y_local[2]],
            [z_local[0], z_local[1], z_local[2]],
        ],
        dtype=float,
    )


def truss_transformation(p0: XYZ, p1: XYZ) -> Matrix6:
    """6x6 transformation for truss translations: u_local = T @ u_global."""
    rot = local_rotation_matrix(p0, p1)
    t = np.zeros((6, 6), dtype=float)
    t[0:3, 0:3] = rot
    t[3:6, 3:6] = rot
    return t


def frame_transformation(p0: XYZ, p1: XYZ) -> Matrix12:
    """12x12 transformation for 3D frame: u_local = T @ u_global."""
    rot = local_rotation_matrix(p0, p1)
    t = np.zeros((12, 12), dtype=float)
    # node i
    t[0:3, 0:3] = rot
    t[3:6, 3:6] = rot
    # node j
    t[6:9, 6:9] = rot
    t[9:12, 9:12] = rot
    return t


def transform_stiffness(local_k: np.ndarray, t: np.ndarray) -> np.ndarray:
    """k_global = T^T * k_local * T"""
    return t.T @ local_k @ t


def transform_load(local_f: np.ndarray, t: np.ndarray) -> np.ndarray:
    """f_global = T^T @ f_local"""
    return t.T @ local_f


def _cross(a: XYZ, b: XYZ) -> XYZ:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(v: XYZ) -> float:
    return float((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5)


def _normalize(v: XYZ) -> XYZ:
    length = _norm(v)
    if length <= 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / length, v[1] / length, v[2] / length)

