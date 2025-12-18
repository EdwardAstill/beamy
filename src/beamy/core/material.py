from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

@dataclass
class Material:
    name: str
    E: float  # Young's modulus
    G: float  # Shear modulus
    Fy: Optional[float] = None  # Yield stress (needed for AISC checks)
    transparency: bool = False  # Used for plotting (affects alpha)

