"""
Lightweight smoke checks for the frame-first Beamy API.

Run directly:
    python testing/analysis.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from sectiony.library import rhs

# Allow running from repo root without installation.
SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from beamy import AnalysisResult, LoadCase, Material, Member, MemberPointForce, Node, Result
from beamy.frame import Frame


def test_result_class() -> None:
    x = np.array([0.0, 1.0, 2.0, 3.0])
    values = np.array([10.0, 20.0, 15.0, 5.0])
    r = Result(x, values)

    assert list(r) == [(0.0, 10.0), (1.0, 20.0), (2.0, 15.0), (3.0, 5.0)]
    assert r.max == 20.0
    assert r.min == 5.0
    assert abs(r.mean - 12.5) < 1e-12
    assert r.range == 15.0
    assert r[0] == (0.0, 10.0)
    assert abs(r.at(0.5) - 15.0) < 1e-12


def test_analysis_result() -> None:
    x = np.array([0.0, 1.0, 2.0])
    action = Result(x, np.array([10.0, 20.0, 15.0]))
    stress = Result(x, np.array([100.0, 200.0, 150.0]))
    displacement = Result(x, np.array([0.0, 0.1, 0.2]))
    ar = AnalysisResult(action, stress, displacement)
    assert ar.action is action
    assert ar.stress is stress
    assert ar.displacement is displacement


def test_frame_analyze_smoke() -> None:
    mat = Material(name="Steel", E=200e9, G=80e9)
    sec = rhs(b=0.1, h=0.2, t=0.005, r=0.0)

    members = [
        Member("M1",
               start=np.array([0.0, 0.0, 0.0]),
               end=np.array([3.0, 0.0, 0.0]),
               section=sec,
               material=mat,
               orientation=np.array([0.0, 1.0, 0.0])),
    ]
    frame = Frame.from_members(members)
    
    # Apply fixed support at start (auto-generated as N0)
    frame.nodes["N0"].support = "111111"

    lc = LoadCase("tip load")
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=3.0,
            force=np.array([0.0, 0.0, -1000.0]),
            coords="global",
            position_type="absolute",
        )
    )

    res = frame.analyze(lc)
    assert res.converged
    assert "N0" in res.nodal_displacements


def main() -> None:
    test_result_class()
    test_analysis_result()
    test_frame_analyze_smoke()
    print("OK")


if __name__ == "__main__":
    main()


