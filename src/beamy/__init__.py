"""
A lightweight 1D beam analysis package.

Provides:
- Beam definitions (material, section, supports)
- Point + distributed loads in local beam coordinates
- Closed-form static + EB bending analysis for a single beam
"""

from .beam import (
    Beam1D,
    Material,
    Section,
    validate_node,
)

from .loads import (
    PointLoad,
    DistributedLoad,
    LoadCase,
)

from .analysis import (
    BeamAnalysisResult,
    analyze_beam_simple_point_load,
)

from . import utils

__all__ = [
    # beam.py
    "Beam1D",
    "Material",
    "Section",
    "validate_node",

    # loads.py
    "PointLoad",
    "DistributedLoad",
    "LoadCase",

    # analysis.py
    "BeamAnalysisResult",
    "analyze_beam_simple_point_load",

    # utils
    "utils",
]
