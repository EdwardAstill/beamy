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

    def add_point_force(self, pf: PointForce):
        self.point_forces.append(pf)

    def add_moment(self, m: Moment):
        self.moments.append(m)

    def add_distributed_force(self, df: DistributedForce):
        self.dist_forces.append(df)

    def convert_dforces_to_pforces(self, n_points: int = 11) -> None:
        """
        Approximate each distributed force by `n_points` equivalent point forces.
        After this, dist_forces is cleared and only point_forces are used.
        """
        if n_points < 2:
            raise ValueError("n_points must be at least 2")

        new_point_forces: list[PointForce] = []

        for df in self.dist_forces:
            x_start = df.start_position
            x_end = df.end_position
            w_start = df.start_force  # [Fx, Fy, Fz] per unit length
            w_end = df.end_force

            ts = np.linspace(0.0, 1.0, n_points)
            pos_interp = np.outer(1 - ts, x_start) + np.outer(ts, x_end)
            w_interp = np.outer(1 - ts, w_start) + np.outer(ts, w_end)

            L_total = float(x_end[0] - x_start[0])
            if L_total == 0.0:
                F_vec = w_start * 0.0
                new_point_forces.append(PointForce(point=x_start, force=F_vec))
                continue

            dx = L_total / (n_points - 1)

            for p, w in zip(pos_interp, w_interp):
                F_vec = w * dx
                new_point_forces.append(PointForce(point=p, force=F_vec))

        self.point_forces.extend(new_point_forces)
        self.dist_forces.clear()

    # --------------------
    # Force resultants
    # --------------------
    @property
    def Fxs(self) -> list[tuple[float, float]]:
        """List of (x, Fx) from point forces only (d-forces must be pre-discretised)."""
        forces: list[tuple[float, float]] = []
        for pf in self.point_forces:
            x = float(pf.point[0])
            Fx = float(pf.force[0])
            if Fx != 0.0:
                forces.append((x, Fx))
        return sorted(forces, key=lambda p: p[0])

    @property
    def Fys(self) -> list[tuple[float, float]]:
        """List of (x, Fy) from point forces only."""
        forces: list[tuple[float, float]] = []
        for pf in self.point_forces:
            x = float(pf.point[0])
            Fy = float(pf.force[1])
            if Fy != 0.0:
                forces.append((x, Fy))
        return sorted(forces, key=lambda p: p[0])

    @property
    def Fzs(self) -> list[tuple[float, float]]:
        """List of (x, Fz) from point forces only."""
        forces: list[tuple[float, float]] = []
        for pf in self.point_forces:
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
        """
        List of (x, Mx) including:
        - torsion from eccentric transverse point forces: Mx = y*Fz - z*Fy
        - explicit torsional moments Moment.moment[0]
        """
        torsions: list[tuple[float, float]] = []

        # From eccentric point forces
        for pf in self.point_forces:
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
        """
        List of (x, My) including:
        - bending from eccentric axial forces: My = z * Fx
        - explicit bending moments about y: Moment.moment[1]
        """
        moments: list[tuple[float, float]] = []

        # From eccentric axial point forces
        for pf in self.point_forces:
            x = float(pf.point[0])
            z = float(pf.point[2])
            Fx = float(pf.force[0])
            My = z * Fx
            if My != 0.0:
                moments.append((x, My))

        # From explicit My moments
        for m in self.moments:
            x = float(m.x)
            My = float(m.moment[1])
            if My != 0.0:
                moments.append((x, My))

        return sorted(moments, key=lambda p: p[0])

    @property
    def Mzs(self) -> list[tuple[float, float]]:
        """
        List of (x, Mz) including:
        - bending from eccentric axial forces: Mz = -y * Fx
        - explicit bending moments about z: Moment.moment[2]
        """
        moments: list[tuple[float, float]] = []

        # From eccentric axial point forces
        for pf in self.point_forces:
            x = float(pf.point[0])
            y = float(pf.point[1])
            Fx = float(pf.force[0])
            Mz = -y * Fx
            if Mz != 0.0:
                moments.append((x, Mz))

        # From explicit Mz moments
        for m in self.moments:
            x = float(m.x)
            Mz = float(m.moment[2])
            if Mz != 0.0:
                moments.append((x, Mz))

        return sorted(moments, key=lambda p: p[0])
