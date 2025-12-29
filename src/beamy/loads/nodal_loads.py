from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Tuple, Union

from beamy.model.node import NodeId

CoordSys = Literal["global", "local"]
ForceVector = Tuple[float, float, float]
MomentVector = Tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class NodalLoad:
    node_id: NodeId
    forces: ForceVector = (0.0, 0.0, 0.0)
    moments: MomentVector = (0.0, 0.0, 0.0)
    coord_sys: CoordSys = "global"

    def __post_init__(self) -> None:
        self._validate_vector(self.forces, "forces")
        self._validate_vector(self.moments, "moments")

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.node_id not in frame.nodes:
            raise ValueError(f"node {self.node_id} not found for nodal load")

    @staticmethod
    def _validate_vector(vector: Tuple[float, float, float], label: str) -> None:
        if len(vector) != 3:
            raise ValueError(f"{label} must have 3 components")
        _ = tuple(float(value) for value in vector)

