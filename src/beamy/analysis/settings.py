from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

SolverKind = Literal["linear", "second_order"]
DimensionKind = Literal["3d"]
FormulationKind = Literal["linear", "pdelta", "corotational"]


@dataclass
class AnalysisSettings:
    dimension: DimensionKind = "3d"
    solver: SolverKind = "linear"
    formulation: FormulationKind = "linear"
    verbose: bool = False
    max_iter: int = 50
    tol: float = 1e-6
    relaxation: float = 1.0

    def validate(self) -> None:
        if self.dimension != "3d":
            raise ValueError("only 3d is supported in the MVP")
        if self.solver not in ("linear", "second_order"):
            raise ValueError("solver must be 'linear' or 'second_order'")
        if self.formulation not in ("linear", "pdelta", "corotational"):
            raise ValueError("formulation must be 'linear', 'pdelta', or 'corotational'")
        if self.max_iter <= 0:
            raise ValueError("max_iter must be positive")
        if self.tol <= 0.0:
            raise ValueError("tol must be positive")
        if self.relaxation <= 0.0 or self.relaxation > 1.0:
            raise ValueError("relaxation must be in (0, 1]")

