from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from beamy.loads.member_loads import MemberLoad
from beamy.loads.nodal_loads import NodalLoad
from beamy.loads.supports import PrescribedDisplacement, Support


@dataclass
class LoadCase:
    """
    A single analysis load case.

    Pure data + validation. No conversion to global vectors here.
    """

    name: str
    nodal_loads: List[NodalLoad] = field(default_factory=list)
    member_loads: List[MemberLoad] = field(default_factory=list)
    supports: List[Support] = field(default_factory=list)
    prescribed_displacements: List[PrescribedDisplacement] = field(default_factory=list)

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        for load in self.nodal_loads:
            load.validate(frame)
        for mload in self.member_loads:
            mload.validate(frame)
        for support in self.supports:
            support.validate(frame)
        for displacement in self.prescribed_displacements:
            displacement.validate(frame)

