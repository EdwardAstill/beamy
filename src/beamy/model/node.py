from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, Tuple, Union, cast

Number = Union[int, float]
XYZ = Tuple[float, float, float]
NodeId = Union[int, str]


@dataclass(slots=True)
class Node:
    """
    Node data structure for 3D frames/trusses.

    Only holds geometry + optional metadata. Validation ensures 3D coords.
    """

    id: NodeId
    xyz: XYZ
    meta: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        coords = self._coerce_xyz(self.xyz)
        object.__setattr__(self, "xyz", coords)

    @staticmethod
    def _coerce_xyz(xyz: Iterable[Number]) -> XYZ:
        coords = tuple(float(value) for value in xyz)
        if len(coords) != 3:
            raise ValueError("xyz must contain exactly 3 coordinates")
        if any(value != value for value in coords):  # NaN check
            raise ValueError("xyz coordinates must be real numbers")
        return cast(XYZ, coords)

