from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from sectiony.library import rhs

# Allow running directly from repo root without installation.
SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from beamy import LoadCase, LoadedMember, Material, MemberPointForce


def _reaction(result, node_id: str) -> np.ndarray:
    if node_id in result.reactions:
        return result.reactions[node_id]
    return np.zeros(6, dtype=float)


def test_cantilever_end_load() -> None:
    """Cantilever with end point load (check internal force magnitudes)."""
    L = 5.0
    P = -1000.0  # N downward in global -Z

    steel = Material(name="Steel", E=200e9, G=80e9)
    section = rhs(b=0.1, h=0.2, t=0.005, r=0.0)

    loads = LoadCase(name="End Load")
    loads.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=L,
            force=np.array([0.0, 0.0, P]),
            coords="global",
            position_type="absolute",
        )
    )

    lm = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111111",  # fixed
        support_end="000000",  # free
        load_case=loads,
    )

    profile = lm.member_demand().actions(points=401)

    # For a tip load, expected magnitudes:
    expected_m0 = abs(P) * L
    expected_v = abs(P)

    m0 = float(profile.bending_y.at(0.0))
    v0 = float(profile.shear_z.at(0.0))
    mL = float(profile.bending_y.at(L))

    assert abs(abs(m0) - expected_m0) / expected_m0 < 0.03
    assert abs(abs(v0) - expected_v) / expected_v < 0.03
    assert abs(mL) < 1e-6 * expected_m0


def test_simply_supported_center_load() -> None:
    """Simply supported with center point load (check moment magnitude + reactions)."""
    L = 10.0
    P = -5000.0  # N downward in global -Z

    steel = Material(name="Steel", E=200e9, G=80e9)
    section = rhs(b=0.15, h=0.25, t=0.006, r=0.0)

    loads = LoadCase(name="Center Load")
    loads.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=L / 2.0,
            force=np.array([0.0, 0.0, P]),
            coords="global",
            position_type="absolute",
        )
    )

    lm = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111100",
        support_end="011000",
        load_case=loads,
    )

    profile = lm.member_demand().actions(points=801)
    m_mid = float(profile.bending_y.at(L / 2.0))

    expected_m_mid = abs(P) * L / 4.0
    assert abs(abs(m_mid) - expected_m_mid) / expected_m_mid < 0.05

    # Reactions (approx): |R1| ~= |R2| ~= |P|/2
    result = lm.analysis_result
    r0 = _reaction(result, "N0")
    r1 = _reaction(result, "N1")
    rz0 = float(r0[2])
    rz1 = float(r1[2])

    expected_r = abs(P) / 2.0
    assert abs(abs(rz0) - expected_r) / expected_r < 0.05
    assert abs(abs(rz1) - expected_r) / expected_r < 0.05


