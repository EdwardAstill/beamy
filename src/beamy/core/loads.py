from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Literal
import numpy as np
from .support import validate_support_type


@dataclass
class NodalForce:
    """Force applied directly at a node in a 3D frame."""
    node_id: str
    force: np.ndarray
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None
    
    def __post_init__(self) -> None:
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)
        if self.force.shape != (3,):
            raise ValueError(f"NodalForce must be a 3D vector [FX, FY, FZ], got shape {self.force.shape}")
        if self.coords not in ("global", "local"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.coords == "local" and not self.reference_member_id:
            raise ValueError("reference_member_id is required when coords='local'")


@dataclass
class NodalMoment:
    """Moment applied directly at a node in a 3D frame."""
    node_id: str
    moment: np.ndarray
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None
    
    def __post_init__(self) -> None:
        if not isinstance(self.moment, np.ndarray):
            self.moment = np.array(self.moment, dtype=float)
        if self.moment.shape != (3,):
            raise ValueError(f"NodalMoment must be a 3D vector [MX, MY, MZ], got shape {self.moment.shape}")
        if self.coords not in ("global", "local"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.coords == "local" and not self.reference_member_id:
            raise ValueError("reference_member_id is required when coords='local'")


@dataclass
class NodalSpring:
    """Elastic nodal spring stiffness attached to ground at a node.

    K is a 6x6 stiffness matrix in the node DOF basis:
        [Ux, Uy, Uz, Rx, Ry, Rz]

    Coordinates:
        - coords='global': K is expressed in global axes
        - coords='local' : K is expressed in the local axes of reference_member_id
          and will be transformed to global during analysis.
    """

    node_id: str
    K: np.ndarray
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None

    def __post_init__(self) -> None:
        if not isinstance(self.K, np.ndarray):
            self.K = np.array(self.K, dtype=float)
        if self.K.shape != (6, 6):
            raise ValueError(f"NodalSpring.K must be shape (6,6), got {self.K.shape}")
        if self.coords not in ("global", "local"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.coords == "local" and not self.reference_member_id:
            raise ValueError("reference_member_id is required when coords='local'")


@dataclass
class MemberPointForce:
    """Point force applied along a member in a 3D frame."""
    member_id: str
    position: float
    force: np.ndarray
    coords: Literal["local", "global"] = "local"
    position_type: Literal["absolute", "relative"] = "absolute"
    
    def __post_init__(self) -> None:
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)
        if self.force.shape != (3,):
            raise ValueError(f"MemberPointForce must be a 3D vector, got shape {self.force.shape}")
        if self.coords not in ("local", "global"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.position_type not in ("absolute", "relative"):
            raise ValueError(f"position_type must be 'absolute' or 'relative', got '{self.position_type}'")
        if self.position_type == "relative" and not (0 <= self.position <= 1):
            raise ValueError(f"Relative position must be in range [0, 1], got {self.position}")
        if self.position < 0:
            raise ValueError(f"Position must be non-negative, got {self.position}")


@dataclass
class MemberPointMoment:
    """Point moment applied along a member in a 3D frame."""
    member_id: str
    position: float
    moment: np.ndarray
    coords: Literal["local", "global"] = "local"
    position_type: Literal["absolute", "relative"] = "absolute"

    def __post_init__(self) -> None:
        if not isinstance(self.moment, np.ndarray):
            self.moment = np.array(self.moment, dtype=float)
        if self.moment.shape != (3,):
            raise ValueError(f"MemberPointMoment must be a 3D vector, got shape {self.moment.shape}")
        if self.coords not in ("local", "global"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.position_type not in ("absolute", "relative"):
            raise ValueError(f"position_type must be 'absolute' or 'relative', got '{self.position_type}'")
        if self.position_type == "relative" and not (0 <= self.position <= 1):
            raise ValueError(f"Relative position must be in range [0, 1], got {self.position}")
        if self.position < 0:
            raise ValueError(f"Position must be non-negative, got {self.position}")


@dataclass
class MemberDistributedForce:
    """Distributed force along a member segment in a 3D frame."""
    member_id: str
    start_position: float
    end_position: Optional[float]  # None means full member length
    start_force: np.ndarray
    end_force: np.ndarray
    coords: Literal["local", "global"] = "local"
    
    def __post_init__(self) -> None:
        if not isinstance(self.start_force, np.ndarray):
            self.start_force = np.array(self.start_force, dtype=float)
        if self.start_force.shape != (3,):
            raise ValueError(f"start_force must be a 3D vector, got shape {self.start_force.shape}")
        if not isinstance(self.end_force, np.ndarray):
            self.end_force = np.array(self.end_force, dtype=float)
        if self.end_force.shape != (3,):
            raise ValueError(f"end_force must be a 3D vector, got shape {self.end_force.shape}")
        if self.coords not in ("local", "global"):
            raise ValueError(f"coords must be 'local' or 'global', got '{self.coords}'")
        if self.start_position < 0:
            raise ValueError(f"start_position must be non-negative, got {self.start_position}")
        if self.end_position is not None and float(self.end_position) < float(self.start_position):
            raise ValueError(
                f"end_position ({self.end_position}) must be >= start_position ({self.start_position}) or None (full length)"
            )


@dataclass(frozen=True)
class MemberPointSupport:
    """A support (restraint to ground) applied at a point along a member.
    
    This is part of the load case boundary conditions, not member geometry.
    Internally, the analysis inserts a node at this position and applies the support.
    
    Support string format is 6 digits for [UX, UY, UZ, RX, RY, RZ].
    """
    member_id: str
    position: float
    support: str
    position_type: Literal["absolute", "relative"] = "absolute"

    def __post_init__(self) -> None:
        validate_support_type(self.support)
        if self.position_type not in ("absolute", "relative"):
            raise ValueError(f"position_type must be 'absolute' or 'relative', got '{self.position_type}'")
        if self.position_type == "relative" and not (0.0 <= float(self.position) <= 1.0):
            raise ValueError(f"Relative position must be in range [0, 1], got {self.position}")
        if float(self.position) < 0.0:
            raise ValueError(f"Position must be non-negative, got {self.position}")


@dataclass(frozen=True)
class MemberSupport:
    """A support applied to all nodes along a member (start, end, and any intermediate nodes).
    
    This is a convenience for applying the same support condition to an entire member.
    """
    member_id: str
    support: str

    def __post_init__(self) -> None:
        validate_support_type(self.support)


@dataclass
class LoadCase:
    """Loads and boundary conditions applied to a structure (standalone member or multi-member frame)."""
    name: str
    nodal_forces: List[NodalForce] = field(default_factory=list)
    nodal_moments: List[NodalMoment] = field(default_factory=list)
    member_point_forces: List[MemberPointForce] = field(default_factory=list)
    member_point_moments: List[MemberPointMoment] = field(default_factory=list)
    member_distributed_forces: List[MemberDistributedForce] = field(default_factory=list)
    nodal_springs: List[NodalSpring] = field(default_factory=list)
    member_point_supports: List[MemberPointSupport] = field(default_factory=list)
    member_supports: List[MemberSupport] = field(default_factory=list)
    
    def add_nodal_force(
        self,
        node_id: str,
        force: np.ndarray,
        coords: str = "global",
        reference_member_id: Optional[str] = None,
    ) -> None:
        self.nodal_forces.append(
            NodalForce(node_id=node_id, force=force, coords=coords, reference_member_id=reference_member_id)
        )
    
    def add_nodal_moment(
        self,
        node_id: str,
        moment: np.ndarray,
        coords: str = "global",
        reference_member_id: Optional[str] = None,
    ) -> None:
        self.nodal_moments.append(
            NodalMoment(node_id=node_id, moment=moment, coords=coords, reference_member_id=reference_member_id)
        )
    
    def add_member_point_force(self, member_id: str, position: float, force: np.ndarray, 
                               coords: str = "local", position_type: str = "absolute") -> None:
        self.member_point_forces.append(MemberPointForce(member_id, position, force, coords, position_type))

    def add_member_point_moment(self, member_id: str, position: float, moment: np.ndarray,
                                coords: str = "local", position_type: str = "absolute") -> None:
        self.member_point_moments.append(MemberPointMoment(member_id, position, moment, coords, position_type))

    def add_member_point_load(
        self,
        member_id: str,
        position: float,
        force: Optional[np.ndarray] = None,
        moment: Optional[np.ndarray] = None,
        coords: str = "local",
        position_type: str = "absolute",
    ) -> None:
        if force is not None:
            self.add_member_point_force(member_id, position, force, coords=coords, position_type=position_type)
        if moment is not None:
            self.add_member_point_moment(member_id, position, moment, coords=coords, position_type=position_type)
    
    def add_member_distributed_force(
        self,
        member_id: str,
        start_position: float,
        end_position: Optional[float],
        start_force: np.ndarray,
        end_force: np.ndarray,
        coords: str = "local",
    ) -> None:
        self.member_distributed_forces.append(
            MemberDistributedForce(member_id, start_position, end_position, start_force, end_force, coords)
        )
    
    def add_member_uniform_force(self, member_id: str, force: np.ndarray, coords: str = "local") -> None:
        self.member_distributed_forces.append(MemberDistributedForce(member_id, 0.0, None, force, force, coords))

    def add_nodal_spring(
        self,
        node_id: str,
        K: np.ndarray,
        coords: str = "global",
        reference_member_id: Optional[str] = None,
    ) -> None:
        self.nodal_springs.append(NodalSpring(node_id=node_id, K=K, coords=coords, reference_member_id=reference_member_id))

    def add_member_point_support(
        self,
        member_id: str,
        position: float,
        support: str,
        position_type: str = "absolute",
    ) -> None:
        """Add a point support (restraint to ground) along a member."""
        self.member_point_supports.append(
            MemberPointSupport(member_id=member_id, position=position, support=support, position_type=position_type)
        )

    def add_member_support(self, member_id: str, support: str) -> None:
        """Add a support to all nodes along a member."""
        self.member_supports.append(MemberSupport(member_id=member_id, support=support))

    def __repr__(self) -> str:
        n_nodal = len(self.nodal_forces) + len(self.nodal_moments) + len(self.nodal_springs)
        n_member = (
            len(self.member_point_forces)
            + len(self.member_point_moments)
            + len(self.member_distributed_forces)
        )
        n_supports = len(self.member_point_supports) + len(self.member_supports)
        return f"LoadCase('{self.name}', {n_nodal} nodal loads, {n_member} member loads, {n_supports} member supports)"
