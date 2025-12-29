from __future__ import annotations

from typing import Tuple

import numpy as np

from beamy.analysis.elements.truss3d import global_stiffness_12x12

XYZ = Tuple[float, float, float]


def cable_stiffness_12x12(p0: XYZ, p1: XYZ, e_modulus: float, area: float) -> np.ndarray:
    """
    Cable stiffness is the same as a truss stiffness when the cable is active.

    Tension-only activation/deactivation is handled by the solver (active-set iteration).
    """
    return global_stiffness_12x12(p0, p1, e_modulus, area)


def cable_is_active(axial_force_value: float) -> bool:
    """Return True if a cable should be active (tension-only)."""
    return float(axial_force_value) >= 0.0

