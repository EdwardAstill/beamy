from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Dict, Iterable, Literal, Tuple, Union, cast

from beamy.model.dof import ReleaseMask6
from beamy.model.node import NodeId, XYZ

MemberKind = Literal["truss", "beam", "cable"]

DEFAULT_RELEASE: ReleaseMask6 = (False, False, False, False, False, False)


@dataclass(slots=True)
class Member:
    """
    Member data structure (3D).

    Supports:
    - truss: axial-only (translations only)
    - beam: 3D Euler-Bernoulli space frame (6 dof/node)
    - cable: tension-only (iterative active-set)
    """

    id: Union[int, str]
    end_i: NodeId
    end_j: NodeId
    kind: MemberKind = "truss"
    area: float = 0.0
    youngs_modulus: float = 0.0
    shear_modulus: float = 0.0
    iy: float = 0.0
    iz: float = 0.0
    j: float = 0.0
    density: float = 0.0
    stations: Tuple[float, ...] = field(default_factory=tuple)
    release_i: ReleaseMask6 = DEFAULT_RELEASE
    release_j: ReleaseMask6 = DEFAULT_RELEASE
    meta: Dict[str, Union[str, float]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._validate_kind()
        self._validate_section()
        self._validate_stations()
        self._validate_releases(self.release_i, "release_i")
        self._validate_releases(self.release_j, "release_j")

    def length(self, frame: "Frame") -> float:
        from beamy.model.frame import Frame  # local import to avoid cycle

        if not isinstance(frame, Frame):
            raise TypeError("frame must be a Frame")
        node_i = self._require_node(frame, self.end_i)
        node_j = self._require_node(frame, self.end_j)
        return self._distance(node_i.xyz, node_j.xyz)

    def local_axis(self, frame: "Frame") -> XYZ:
        from beamy.analysis.transformations import unit_vector

        node_i = self._require_node(frame, self.end_i)
        node_j = self._require_node(frame, self.end_j)
        return unit_vector(node_i.xyz, node_j.xyz)

    def _validate_kind(self) -> None:
        if self.kind not in ("truss", "beam", "cable"):
            raise ValueError("kind must be 'truss', 'beam', or 'cable'")

    def _validate_section(self) -> None:
        if self.area <= 0.0:
            raise ValueError("area must be positive")
        if self.youngs_modulus <= 0.0:
            raise ValueError("youngs_modulus must be positive")
        if self.kind == "beam":
            if self.shear_modulus <= 0.0:
                raise ValueError("shear_modulus must be positive for beam members")
            if self.iy <= 0.0 or self.iz <= 0.0 or self.j <= 0.0:
                raise ValueError("iy, iz, and j must be positive for beam members")
        if self.density < 0.0:
            raise ValueError("density cannot be negative")

    def _validate_stations(self) -> None:
        if not self.stations:
            return
        ordered = sorted(self.stations)
        if any(value < 0.0 or value > 1.0 for value in ordered):
            raise ValueError("stations must lie in [0, 1]")
        object.__setattr__(self, "stations", tuple(ordered))

    @staticmethod
    def _validate_releases(releases: ReleaseMask6, label: str) -> None:
        if len(releases) != 6:
            raise ValueError(f"{label} must have 6 booleans")

    @staticmethod
    def _distance(p0: XYZ, p1: XYZ) -> float:
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        dz = p1[2] - p0[2]
        return float((dx * dx + dy * dy + dz * dz) ** 0.5)

    @staticmethod
    def _require_node(frame: "Frame", node_id: NodeId) -> "Node":
        if node_id not in frame.nodes:
            raise KeyError(f"node {node_id} not found in frame")
        return cast("Node", frame.nodes[node_id])


if TYPE_CHECKING:
    from beamy.model.frame import Frame
    from beamy.model.node import Node

