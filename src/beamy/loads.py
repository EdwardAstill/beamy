# loads.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List
import numpy as np

@dataclass
class PointForce:
    """
    Point force applied to a 1D beam.

    point : np.ndarray([x, y, z])
        Position in *local* coordinates.
        x = station along beam axis
        y,z = offset from shear centre

    force : np.ndarray([Fx, Fy, Fz])
        Force vector in local coordinates.
    """
    point: np.ndarray          # shape (3,) -> [x, y, z]
    force: np.ndarray          # shape (3,) -> [Fx, Fy, Fz]

    @property
    def x(self) -> float:
        """Station along beam axis (0 <= x <= L)."""
        return float(self.point[0])

    @property
    def offset(self) -> np.ndarray:
        """
        Offset from shear centre at that station.
        (Only y,z matter for eccentricity.)
        """
        return self.point  # [x, y, z], but solver will use [y, z]

    def equivalent_torsion(self) -> float:
        """
        Compute torsion from eccentricity: (r × F)_x.
        r = [0, y, z]
        """
        r = np.array([0.0, self.point[1], self.point[2]])
        M = np.cross(r, self.force)
        return float(M[0])  # torsion about x axis

@dataclass
class Moment:
    """
    Pure moment applied at a station along the beam axis.

    x      : station along beam axis (0 <= x <= L)
    moment : np.ndarray([T, My, Mz]) in local coordinates
    """
    x: float
    moment: np.ndarray  # shape (3,) -> [T, My, Mz]


@dataclass
class DistributedForce:
    """
    Uniform distributed load per unit length over [x_start, x_end] in local axes.
    w: [n_x, v_y, v_z] (force per unit length).
    """
    start_position: np.ndarray
    end_position: np.ndarray
    start_force: np.ndarray
    end_force: np.ndarray #[Fx, Fy, Fz])

@dataclass
class LoadCase:
    name: str
    point_forces: List[PointForce] = field(default_factory=list)
    moments: List[Moment] = field(default_factory=list)
    dist_forces: List[DistributedForce] = field(default_factory=list)

    def add_point_force(self, pf: PointForce):
        self.point_forces.append(pf)

    def add_point_moment(self, pm: Moment):
        self.moments.append(pm)

    def add_distributed_force(self, df: DistributedForce):
        self.dist_forces.append(df)
    
    def get_axial_forces(self, x: float) -> float:
        """
        converts distributed forces to point forces
        returns a list of tuples (x, Fx) for the axial forces at the given x positions
        """
        axial_forces = []
        for pf in self.point_forces:
            if pf.x == x:
                axial_forces.append((x, pf.force[0]))
        for df in self.dist_forces:
            if df.x_start <= x <= df.x_end:
                axial_forces.append((x, df.w[0] * (x - df.x_start)))
        return axial_forces

