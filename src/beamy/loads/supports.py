from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, Union

from beamy.model.dof import RestraintMask6
from beamy.model.node import NodeId

Disp6 = Tuple[float, float, float, float, float, float]
K6 = Tuple[float, float, float, float, float, float]


@dataclass(frozen=True, slots=True)
class Support:
    target: NodeId
    restrained_dofs: RestraintMask6

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.target not in frame.nodes:
            raise ValueError(f"support references missing node {self.target}")
        if len(self.restrained_dofs) != 6:
            raise ValueError("restrained_dofs must have length 6")


@dataclass(frozen=True, slots=True)
class Spring:
    target: NodeId
    k_trans: K6
    k_rot: K6

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.target not in frame.nodes:
            raise ValueError(f"spring references missing node {self.target}")


@dataclass(frozen=True, slots=True)
class PrescribedDisplacement:
    target: NodeId
    values: Disp6

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.target not in frame.nodes:
            raise ValueError(f"prescribed displacement references missing node {self.target}")
        if len(self.values) != 6:
            raise ValueError("values must have 6 components")

