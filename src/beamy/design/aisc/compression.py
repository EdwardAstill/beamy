from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Literal

DesignMethod = Literal["LRFD", "ASD"]


@dataclass(frozen=True)
class CompressionResult:
    method: DesignMethod
    area: float
    fy: float
    e: float
    k_l_over_r: float
    fe: float
    fcr: float
    pn: float
    strength: float
    phi: float
    omega: float
    limit_state: str


def aisc360_e3_fcr(k_l_over_r: float, fy: float, e: float = 29000.0) -> CompressionResult:
    """
    AISC 360-16 (and later) Chapter E, E3 flexural buckling column curve.

    This is a simplified implementation:
    - Uses only the E3 column curve based on KL/r.
    - Does not cover torsional / flexural-torsional buckling (E4) or slender elements (E7).
    """
    if k_l_over_r <= 0.0:
        raise ValueError("k_l_over_r must be positive")
    if fy <= 0.0:
        raise ValueError("fy must be positive")
    if e <= 0.0:
        raise ValueError("e must be positive")

    fe = (pi * pi * e) / (k_l_over_r * k_l_over_r)
    slender_limit = 4.71 * sqrt(e / fy)

    if k_l_over_r <= slender_limit:
        fcr = (0.658 ** (fy / fe)) * fy
        limit_state = "E3 inelastic flexural buckling"
    else:
        fcr = 0.877 * fe
        limit_state = "E3 elastic flexural buckling"

    # area and method-specific strength are handled in a separate function,
    # but we return a structured result here for reuse.
    return CompressionResult(
        method="LRFD",
        area=0.0,
        fy=fy,
        e=e,
        k_l_over_r=k_l_over_r,
        fe=fe,
        fcr=fcr,
        pn=0.0,
        strength=0.0,
        phi=0.9,
        omega=1.67,
        limit_state=limit_state,
    )


def aisc360_compression_capacity(
    area: float,
    k_l_over_r: float,
    fy: float,
    e: float = 29000.0,
    method: DesignMethod = "LRFD",
) -> CompressionResult:
    """
    Compute AISC compression capacity for flexural buckling (E3).

    Returns:
    - pn: nominal compression strength
    - strength: phi*Pn (LRFD) or Pn/Omega (ASD)
    """
    if area <= 0.0:
        raise ValueError("area must be positive")
    curve = aisc360_e3_fcr(k_l_over_r, fy, e)
    pn = curve.fcr * area

    if method == "LRFD":
        phi = 0.9
        omega = 1.67
        strength = phi * pn
    else:
        phi = 0.9
        omega = 1.67
        strength = pn / omega

    return CompressionResult(
        method=method,
        area=area,
        fy=fy,
        e=e,
        k_l_over_r=k_l_over_r,
        fe=curve.fe,
        fcr=curve.fcr,
        pn=pn,
        strength=strength,
        phi=phi,
        omega=omega,
        limit_state=curve.limit_state,
    )

