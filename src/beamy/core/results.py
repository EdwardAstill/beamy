from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class Result:
    """
    Wraps analysis results (x, values) to provide convenient accessors.
    Used for force, stress, and displacement distributions.
    """
    _x: np.ndarray
    _values: np.ndarray

    def __iter__(self):
        return zip(self._x.tolist(), self._values.tolist())

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return list(self)[idx]
        return (self._x[idx], self._values[idx])

    @property
    def max(self) -> float:
        return float(np.max(self._values))

    @property
    def min(self) -> float:
        return float(np.min(self._values))

    @property
    def abs_max(self) -> float:
        return float(np.max(np.abs(self._values)))

    @property
    def mean(self) -> float:
        return float(np.mean(self._values))

    @property
    def range(self) -> float:
        return float(np.ptp(self._values))

    def at(self, x_loc: float) -> float:
        """Interpolate the value at a specific position."""
        return float(np.interp(x_loc, self._x, self._values))

@dataclass
class AnalysisResult:
    """Container for a triple of action, stress, and displacement results."""
    action: Result
    stress: Result
    displacement: Result

