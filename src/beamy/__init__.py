"""
beamy - A lightweight 3D linear frame/truss analysis library (MVP).

Only the easy parts are implemented for now:
- 3D truss/beam-ready data structures
- Linear solver pipeline (mesh → assembly → solve → recovery for truss)
- Simple plotting/export stubs

Hard parts intentionally left out for later:
- AISC design checks
- Second-order analysis
- Cable/tension-only elements
- Full 3D beam formulation
"""

from beamy.analysis import run_analysis
from beamy.analysis.settings import AnalysisSettings
from beamy.loads.loadcase import LoadCase
from beamy.model.frame import Frame
from beamy.model.member import Member
from beamy.model.node import Node
from beamy.results.frame_result import FrameResult
from beamy.results.member_result import MemberResult

__all__ = [
    "Frame",
    "Member",
    "Node",
    "LoadCase",
    "AnalysisSettings",
    "FrameResult",
    "MemberResult",
    "run_analysis",
]

