# utils.py
from __future__ import annotations
import numpy as np

def heaviside(x: float) -> float:
    """Simple Heaviside step (0 for x<0, 1 for x>=0)."""
    return 0.0 if x < 0.0 else 1.0

def clip_to_span(x: float, L: float) -> float:
    """Clamp x to [0, L] for safety."""
    return max(0.0, min(L, x))
