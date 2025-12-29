from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class InteractionResult:
    ratio: float
    governing: str
    axial_ratio: float
    mx_ratio: float
    my_ratio: float

    @property
    def ok(self) -> bool:
        return self.ratio <= 1.0


def aisc360_h1_interaction(
    pu: float,
    pc: float,
    mux: float,
    mcx: float,
    muy: float,
    mcy: float,
) -> InteractionResult:
    """
    AISC 360 Chapter H, H1 interaction (simplified).

    Inputs should be consistent with the design method chosen upstream:
    - pc: available axial strength (phi*Pn or Pn/Omega)
    - mcx/mcy: available flexural strengths (phi*Mn or Mn/Omega)
    """
    if pc <= 0.0:
        raise ValueError("pc must be positive")
    if mcx <= 0.0 or mcy <= 0.0:
        raise ValueError("mcx and mcy must be positive")

    pr_pc = float(pu) / float(pc)
    mx = abs(float(mux)) / float(mcx)
    my = abs(float(muy)) / float(mcy)

    if pr_pc >= 0.2:
        ratio = pr_pc + (8.0 / 9.0) * (mx + my)
        governing = "H1-1a"
    else:
        ratio = (pr_pc / 2.0) + (mx + my)
        governing = "H1-1b"

    return InteractionResult(
        ratio=float(ratio),
        governing=governing,
        axial_ratio=float(pr_pc),
        mx_ratio=float(mx),
        my_ratio=float(my),
    )

