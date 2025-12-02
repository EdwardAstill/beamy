from __future__ import annotations
from dataclasses import dataclass, field
from typing import List
import numpy as np


@dataclass
class PointForce:
    point: np.ndarray   # [x, y, z]
    force: np.ndarray   # [Fx, Fy, Fz]


@dataclass
class Moment:
    x: float
    moment: np.ndarray  # [T, My, Mz]


@dataclass
class DistributedForce:
    start_position: np.ndarray
    end_position: np.ndarray
    start_force: np.ndarray   # [Fx, Fy, Fz] per unit length
    end_force: np.ndarray


@dataclass
class LoadCase:
    """Applied loads only, not considering support reactions."""
    name: str
    point_forces: List[PointForce] = field(default_factory=list)
    moments: List[Moment] = field(default_factory=list)
    dist_forces: List[DistributedForce] = field(default_factory=list)
    dist_force_resolution: int = 11  # Number of points to approximate distributed loads

    def add_point_force(self, pf: PointForce):
        self.point_forces.append(pf)

    def add_moment(self, m: Moment):
        self.moments.append(m)

    def add_distributed_force(self, df: DistributedForce):
        self.dist_forces.append(df)

    @property
    def _effective_point_forces(self) -> List[PointForce]:
        """
        Combines explicit point forces with discretized distributed forces.
        Does NOT modify the internal state of dist_forces.
        """
        # Start with explicit point forces
        all_forces = list(self.point_forces)

        # Discretize distributed forces
        n_points = self.dist_force_resolution
        if n_points < 2:
            n_points = 2

        for df in self.dist_forces:
            x_start = df.start_position
            x_end = df.end_position
            w_start = df.start_force
            w_end = df.end_force

            ts = np.linspace(0.0, 1.0, n_points)
            pos_interp = np.outer(1 - ts, x_start) + np.outer(ts, x_end)
            w_interp = np.outer(1 - ts, w_start) + np.outer(ts, w_end)

            L_total = float(x_end[0] - x_start[0])
            if L_total == 0.0:
                continue

            dx = L_total / (n_points - 1)

            for p, w in zip(pos_interp, w_interp):
                F_vec = w * dx
                all_forces.append(PointForce(point=p, force=F_vec))
        
        return all_forces

    # --------------------
    # Force resultants
    # --------------------
    @property
    def Fxs(self) -> list[tuple[float, float]]:
        """List of (x, Fx) from all forces (explicit + discretized)."""
        forces: list[tuple[float, float]] = []
        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            Fx = float(pf.force[0])
            if Fx != 0.0:
                forces.append((x, Fx))
        return sorted(forces, key=lambda p: p[0])

    @property
    def Fys(self) -> list[tuple[float, float]]:
        """List of (x, Fy) from all forces."""
        forces: list[tuple[float, float]] = []
        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            Fy = float(pf.force[1])
            if Fy != 0.0:
                forces.append((x, Fy))
        return sorted(forces, key=lambda p: p[0])

    @property
    def Fzs(self) -> list[tuple[float, float]]:
        """List of (x, Fz) from all forces."""
        forces: list[tuple[float, float]] = []
        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            Fz = float(pf.force[2])
            if Fz != 0.0:
                forces.append((x, Fz))
        return sorted(forces, key=lambda p: p[0])

    # --------------------
    # Moment resultants
    # --------------------
    @property
    def Mxs(self) -> list[tuple[float, float]]:
        """List of (x, Mx) including torsion from eccentric forces."""
        torsions: list[tuple[float, float]] = []

        # From all point forces (explicit + discretized)
        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            y = float(pf.point[1])
            z = float(pf.point[2])
            Fx, Fy, Fz = map(float, pf.force)
            Mx = y * Fz - z * Fy
            if Mx != 0.0:
                torsions.append((x, Mx))

        # From explicit torsional moments
        for m in self.moments:
            x = float(m.x)
            T = float(m.moment[0])
            if T != 0.0:
                torsions.append((x, T))

        return sorted(torsions, key=lambda p: p[0])

    @property
    def Mys(self) -> list[tuple[float, float]]:
        """List of (x, My) including bending from eccentric forces."""
        moments: list[tuple[float, float]] = []

        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            z = float(pf.point[2])
            Fx = float(pf.force[0])
            My = z * Fx
            if My != 0.0:
                moments.append((x, My))

        for m in self.moments:
            x = float(m.x)
            My = float(m.moment[1])
            if My != 0.0:
                moments.append((x, My))

        return sorted(moments, key=lambda p: p[0])

    @property
    def Mzs(self) -> list[tuple[float, float]]:
        """List of (x, Mz) including bending from eccentric forces."""
        moments: list[tuple[float, float]] = []

        for pf in self._effective_point_forces:
            x = float(pf.point[0])
            y = float(pf.point[1])
            Fx = float(pf.force[0])
            Mz = -y * Fx
            if Mz != 0.0:
                moments.append((x, Mz))

        for m in self.moments:
            x = float(m.x)
            Mz = float(m.moment[2])
            if Mz != 0.0:
                moments.append((x, Mz))

        return sorted(moments, key=lambda p: p[0])
