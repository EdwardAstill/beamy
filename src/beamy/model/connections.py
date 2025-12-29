from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, Tuple, Union

EndLabel = Literal["i", "j"]
ReleaseMask = Tuple[bool, bool, bool, bool, bool, bool]


@dataclass(frozen=True, slots=True)
class EndRef:
    member_id: Union[int, str]
    end: EndLabel


@dataclass(frozen=True, slots=True)
class StationRef:
    member_id: Union[int, str]
    s: float

    def __post_init__(self) -> None:
        if self.s < 0.0 or self.s > 1.0:
            raise ValueError("station s must lie in [0, 1]")


@dataclass(frozen=True, slots=True)
class EndToStation:
    end_ref: EndRef
    station_ref: StationRef
    dofs: Optional[ReleaseMask] = None

    def __post_init__(self) -> None:
        if self.dofs is not None and len(self.dofs) != 6:
            raise ValueError("dofs mask must have length 6")


Connection = Union[EndToStation]

