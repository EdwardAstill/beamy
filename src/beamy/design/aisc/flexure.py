from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

DesignMethod = Literal["LRFD", "ASD"]
Axis = Literal["x", "y"]


@dataclass(frozen=True)
class FlexureResult:
    method: DesignMethod
    axis: Axis
    fy: float
    z: float
    s: float
    mn: float
    strength: float
    phi: float
    omega: float
    limit_state: str


def aisc360_flexure_yielding_capacity(
    fy: float,
    axis: Axis,
    z: Optional[float] = None,
    s: Optional[float] = None,
    method: DesignMethod = "LRFD",
) -> FlexureResult:
    """
    Simplified AISC Chapter F flexure strength (yielding only).

    - If plastic modulus Z is provided, uses Mp = Fy*Z.
    - Else uses elastic section modulus S, My = Fy*S.

    Does NOT include:
    - LTB (lateral-torsional buckling)
    - flange/web local buckling reductions
    - shape-specific Chapter F branching
    """
    if fy <= 0.0:
        raise ValueError("fy must be positive")
    if axis not in ("x", "y"):
        raise ValueError("axis must be 'x' or 'y'")

    z_val = float(z) if z is not None else 0.0
    s_val = float(s) if s is not None else 0.0
    if z_val <= 0.0 and s_val <= 0.0:
        raise ValueError("provide z (plastic modulus) or s (elastic section modulus)")

    if z_val > 0.0:
        mn = fy * z_val
        limit_state = "F1 yielding (plastic)"
    else:
        mn = fy * s_val
        limit_state = "F1 yielding (elastic)"

    phi = 0.9
    omega = 1.67
    strength = phi * mn if method == "LRFD" else mn / omega

    return FlexureResult(
        method=method,
        axis=axis,
        fy=fy,
        z=z_val,
        s=s_val,
        mn=mn,
        strength=strength,
        phi=phi,
        omega=omega,
        limit_state=limit_state,
    )

