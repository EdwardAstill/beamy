from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Tuple, Union

CoordSys = Literal["global", "local"]
Vec3 = Tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class MemberPointLoad:
    member_id: Union[int, str]
    s: float
    vector: Vec3
    coord_sys: CoordSys = "local"

    def __post_init__(self) -> None:
        _validate_station(self.s)
        _validate_vec(self.vector, "vector")

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.member_id not in frame.members:
            raise ValueError(f"member {self.member_id} not found for member point load")


@dataclass(frozen=True, slots=True)
class MemberPointMoment:
    member_id: Union[int, str]
    s: float
    vector: Vec3
    coord_sys: CoordSys = "local"

    def __post_init__(self) -> None:
        _validate_station(self.s)
        _validate_vec(self.vector, "vector")

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.member_id not in frame.members:
            raise ValueError(f"member {self.member_id} not found for member point moment")


@dataclass(frozen=True, slots=True)
class MemberDistributedLoad:
    member_id: Union[int, str]
    s0: float
    s1: float
    w0: Vec3
    w1: Vec3
    coord_sys: CoordSys = "local"

    def __post_init__(self) -> None:
        _validate_station(self.s0)
        _validate_station(self.s1)
        if self.s1 < self.s0:
            raise ValueError("s1 must be >= s0 for distributed load")
        _validate_vec(self.w0, "w0")
        _validate_vec(self.w1, "w1")

    def validate(self, frame: "Frame") -> None:
        from beamy.model.frame import Frame

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        if self.member_id not in frame.members:
            raise ValueError(f"member {self.member_id} not found for member distributed load")


MemberLoad = Union[MemberPointLoad, MemberPointMoment, MemberDistributedLoad]


def _validate_station(value: float) -> None:
    if value < 0.0 or value > 1.0:
        raise ValueError("station must lie in [0, 1]")


def _validate_vec(vector: Vec3, label: str) -> None:
    if len(vector) != 3:
        raise ValueError(f"{label} must have 3 components")
    _ = tuple(float(item) for item in vector)

