from __future__ import annotations
from dataclasses import dataclass
from typing import List
from sectiony import Section

from ..core.material import Material
from ..core.support import Support, validate_support_pairs

@dataclass
class Beam1D:
    """
    Straight prismatic beam along local x in [0, L].
    """
    L: float
    material: Material
    section: Section
    supports: List[Support]

    def __post_init__(self):
        validate_support_pairs(self.supports)

