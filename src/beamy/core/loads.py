from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Any
import numpy as np


@dataclass
class PointForce:
    """A point force applied at a coordinate [x, y, z]. Used primarily for 1D beams."""
    point: np.ndarray   # [x, y, z]
    force: np.ndarray   # [Fx, Fy, Fz]

    def __post_init__(self):
        if not isinstance(self.point, np.ndarray):
            self.point = np.array(self.point, dtype=float)
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)


@dataclass
class Moment:
    """A moment applied at position x along a beam axis."""
    x: float
    moment: np.ndarray  # [T, My, Mz]

    def __post_init__(self):
        if not isinstance(self.moment, np.ndarray):
            self.moment = np.array(self.moment, dtype=float)


@dataclass
class DistributedForce:
    """A distributed force between two points. Used primarily for 1D beams."""
    start_position: np.ndarray
    end_position: np.ndarray
    start_force: np.ndarray   # [Fx, Fy, Fz] per unit length
    end_force: np.ndarray

    def __post_init__(self):
        if not isinstance(self.start_position, np.ndarray):
            self.start_position = np.array(self.start_position, dtype=float)
        if not isinstance(self.end_position, np.ndarray):
            self.end_position = np.array(self.end_position, dtype=float)
        if not isinstance(self.start_force, np.ndarray):
            self.start_force = np.array(self.start_force, dtype=float)
        if not isinstance(self.end_force, np.ndarray):
            self.end_force = np.array(self.end_force, dtype=float)


@dataclass
class NodalForce:
    """Force applied directly at a node in a 3D frame."""
    node_id: str
    force: np.ndarray
    
    def __post_init__(self) -> None:
        if not isinstance(self.force, np.ndarray):
            self.force = np.array(self.force, dtype=float)
        if self.force.shape != (3,):
            raise ValueError(f"NodalForce must be a 3D vector [FX, FY, FZ], got shape {self.force.shape}")


@dataclass
class NodalMoment:
    """Moment applied directly at a node in a 3D frame."""
    node_id: str
    moment: np.ndarray
    
    def __post_init__(self) -> None:
        if not isinstance(self.moment, np.ndarray):
            self.moment = np.array(self.moment, dtype=float)
        if self.moment.shape != (3,):
            raise ValueError(f"NodalMoment must be a 3D vector [MX, MY, MZ], got shape {self.moment.shape}")


@dataclass
class MemberPointForce:
    """Point force applied along a member in a 3D frame."""
    member_id: str
    position: float
    force: np.ndarray
    coords: str = "local"
    position_type: str = "absolute"
    
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
class MemberDistributedForce:
    """Distributed force along a member segment in a 3D frame."""
    member_id: str
    start_position: float
    end_position: float
    start_force: np.ndarray
    end_force: np.ndarray
    coords: str = "local"
    
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
        if self.end_position != -1.0 and self.end_position < self.start_position:
            raise ValueError(f"end_position ({self.end_position}) must be >= start_position ({self.start_position}) or -1 (full length)")


@dataclass
class LoadCase:
    """Applied loads for a 1D beam analysis."""
    name: str
    point_forces: List[PointForce] = field(default_factory=list)
    moments: List[Moment] = field(default_factory=list)
    dist_forces: List[DistributedForce] = field(default_factory=list)
    dist_force_resolution: int = 11

    def add_point_force(self, pf: PointForce):
        self.point_forces.append(pf)

    def add_moment(self, m: Moment):
        self.moments.append(m)

    def add_distributed_force(self, df: DistributedForce):
        self.dist_forces.append(df)

    @property
    def _effective_point_forces(self) -> List[PointForce]:
        all_forces = list(self.point_forces)
        n_points = max(2, self.dist_force_resolution)

        for df in self.dist_forces:
            x_start, x_end = df.start_position, df.end_position
            w_start, w_end = df.start_force, df.end_force

            ts = np.linspace(0.0, 1.0, n_points)
            pos_interp = np.outer(1 - ts, x_start) + np.outer(ts, x_end)
            w_interp = np.outer(1 - ts, w_start) + np.outer(ts, w_end)

            L_total = float(x_end[0] - x_start[0])
            if L_total == 0.0: continue

            dx = L_total / (n_points - 1)
            for i, (p, w) in enumerate(zip(pos_interp, w_interp)):
                weight = 0.5 if (i == 0 or i == n_points - 1) else 1.0
                all_forces.append(PointForce(point=p, force=w * dx * weight))
        
        return all_forces

    @property
    def Fxs(self) -> list[tuple[float, float]]:
        return sorted([(float(pf.point[0]), float(pf.force[0])) 
                      for pf in self._effective_point_forces if pf.force[0] != 0], key=lambda x: x[0])

    @property
    def Fys(self) -> list[tuple[float, float]]:
        return sorted([(float(pf.point[0]), float(pf.force[1])) 
                      for pf in self._effective_point_forces if pf.force[1] != 0], key=lambda x: x[0])

    @property
    def Fzs(self) -> list[tuple[float, float]]:
        return sorted([(float(pf.point[0]), float(pf.force[2])) 
                      for pf in self._effective_point_forces if pf.force[2] != 0], key=lambda x: x[0])

    @property
    def Mxs(self) -> list[tuple[float, float]]:
        torsions = []
        for pf in self._effective_point_forces:
            t = float(pf.point[1] * pf.force[2] - pf.point[2] * pf.force[1])
            torsions.append((float(pf.point[0]), t))
        for m in self.moments:
            torsions.append((float(m.x), float(m.moment[0])))
        return sorted([(x, t) for x, t in torsions if t != 0], key=lambda x: x[0])

    @property
    def Mys(self) -> list[tuple[float, float]]:
        moments = [(float(pf.point[0]), float(pf.point[2] * pf.force[0])) for pf in self._effective_point_forces]
        moments += [(float(m.x), float(m.moment[1])) for m in self.moments]
        return sorted([(x, m) for x, m in moments if m != 0], key=lambda x: x[0])

    @property
    def Mzs(self) -> list[tuple[float, float]]:
        moments = [(float(pf.point[0]), -float(pf.point[1] * pf.force[0])) for pf in self._effective_point_forces]
        moments += [(float(m.x), float(m.moment[2])) for m in self.moments]
        return sorted([(x, m) for x, m in moments if m != 0], key=lambda x: x[0])


@dataclass
class FrameLoadCase:
    """Loads applied to a 3D frame structure."""
    name: str
    nodal_forces: List[NodalForce] = field(default_factory=list)
    nodal_moments: List[NodalMoment] = field(default_factory=list)
    member_point_forces: List[MemberPointForce] = field(default_factory=list)
    member_distributed_forces: List[MemberDistributedForce] = field(default_factory=list)
    
    def add_nodal_force(self, node_id: str, force: np.ndarray) -> None:
        self.nodal_forces.append(NodalForce(node_id=node_id, force=force))
    
    def add_nodal_moment(self, node_id: str, moment: np.ndarray) -> None:
        self.nodal_moments.append(NodalMoment(node_id=node_id, moment=moment))
    
    def add_member_point_force(self, member_id: str, position: float, force: np.ndarray, 
                               coords: str = "local", position_type: str = "absolute") -> None:
        self.member_point_forces.append(MemberPointForce(member_id, position, force, coords, position_type))
    
    def add_member_distributed_force(self, member_id: str, start_position: float, end_position: float, 
                                     start_force: np.ndarray, end_force: np.ndarray, coords: str = "local") -> None:
        self.member_distributed_forces.append(MemberDistributedForce(member_id, start_position, end_position, start_force, end_force, coords))
    
    def add_member_uniform_force(self, member_id: str, force: np.ndarray, coords: str = "local") -> None:
        self.member_distributed_forces.append(MemberDistributedForce(member_id, 0.0, -1.0, force, force, coords))

    def __repr__(self) -> str:
        n_nodal = len(self.nodal_forces) + len(self.nodal_moments)
        n_member = len(self.member_point_forces) + len(self.member_distributed_forces)
        return f"FrameLoadCase('{self.name}', {n_nodal} nodal loads, {n_member} member loads)"
