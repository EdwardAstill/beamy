"""
Standalone member analysis (vNext).

`LoadedMember` is the first-class “beam” object: a single member analyzed in isolation.
It is implemented by building a 2-node `Frame` with one `Member` and delegating to the
frame solver (`Frame.analyze(...)`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, List, Any

import numpy as np
from sectiony import Section

from ..core.material import Material
from ..core.support import validate_support_type
from ..core.loads import LoadCase, MemberPointSupport
from ..core.results import Result, AnalysisResult
from ..frame.node import Node
from ..frame.member import Member
from ..frame.frame import Frame
from ..frame.analysis import FrameAnalysisSettings, FrameAnalysisResult
from ..frame.results import MemberDemand


Axis = Literal["y", "z"]


@dataclass
class LoadedMember:
    """Analyze a single member in isolation under a `LoadCase`.

    Args:
        id: Member identifier (used by member loads in `load_case`)
        start: Global start point (x,y,z)
        end: Global end point (x,y,z)
        section: Cross-section properties (sectiony)
        material: Material properties
        orientation: Local Y axis direction (3-vector)
        support_start: 6-digit support code at start node [UX, UY, UZ, RX, RY, RZ]
        support_end: 6-digit support code at end node [UX, UY, UZ, RX, RY, RZ]
        load_case: Applied loads (frame-style LoadCase). Use load_case.add_member_point_support() for intermediate supports.
        settings: Optional analysis settings
        element_type: "beam" | "truss" | "cable"
        releases: Optional 12-digit release string (member end releases)
        constraints: Optional 12-digit constraint string (member end constraints)
    """

    id: str
    start: np.ndarray
    end: np.ndarray
    section: Section
    material: Material
    orientation: np.ndarray
    support_start: str
    support_end: str
    load_case: LoadCase

    settings: Optional[FrameAnalysisSettings] = None
    element_type: Literal["beam", "truss", "cable"] = "beam"
    releases: Optional[str] = None
    constraints: Optional[str] = None

    _frame: Frame = field(init=False, repr=False)
    _start_node_id: str = field(init=False, default="N0", repr=False)
    _end_node_id: str = field(init=False, default="N1", repr=False)

    def __post_init__(self) -> None:
        self.start = np.asarray(self.start, dtype=float)
        self.end = np.asarray(self.end, dtype=float)
        self.orientation = np.asarray(self.orientation, dtype=float)

        if self.start.shape != (3,):
            raise ValueError(f"start must be a 3D vector, got shape {self.start.shape}")
        if self.end.shape != (3,):
            raise ValueError(f"end must be a 3D vector, got shape {self.end.shape}")
        if self.orientation.shape != (3,):
            raise ValueError(f"orientation must be a 3D vector, got shape {self.orientation.shape}")

        self.support_start = validate_support_type(self.support_start)
        self.support_end = validate_support_type(self.support_end)

        self._validate_load_case()
        self._build_and_analyze()

    def _validate_load_case(self) -> None:
        # Member-targeted loads must reference this member.
        for mpf in self.load_case.member_point_forces:
            if mpf.member_id != self.id:
                raise ValueError(
                    f"Standalone LoadedMember '{self.id}' cannot accept MemberPointForce for member '{mpf.member_id}'"
                )
        for mpm in getattr(self.load_case, "member_point_moments", []):
            if mpm.member_id != self.id:
                raise ValueError(
                    f"Standalone LoadedMember '{self.id}' cannot accept MemberPointMoment for member '{mpm.member_id}'"
                )
        for mdf in self.load_case.member_distributed_forces:
            if mdf.member_id != self.id:
                raise ValueError(
                    f"Standalone LoadedMember '{self.id}' cannot accept MemberDistributedForce for member '{mdf.member_id}'"
                )
        for mps in getattr(self.load_case, "member_point_supports", []):
            if mps.member_id != self.id:
                raise ValueError(
                    f"Standalone LoadedMember '{self.id}' cannot accept MemberPointSupport for member '{mps.member_id}'"
                )

        # Nodal loads are allowed only at the member end nodes.
        for nf in self.load_case.nodal_forces:
            if nf.node_id not in (self._start_node_id, self._end_node_id):
                raise ValueError(
                    f"Standalone LoadedMember nodal force node_id must be '{self._start_node_id}' or '{self._end_node_id}', got '{nf.node_id}'"
                )
        for nm in self.load_case.nodal_moments:
            if nm.node_id not in (self._start_node_id, self._end_node_id):
                raise ValueError(
                    f"Standalone LoadedMember nodal moment node_id must be '{self._start_node_id}' or '{self._end_node_id}', got '{nm.node_id}'"
                )
        for ns in getattr(self.load_case, "nodal_springs", []):
            if ns.node_id not in (self._start_node_id, self._end_node_id):
                raise ValueError(
                    f"Standalone LoadedMember nodal spring node_id must be '{self._start_node_id}' or '{self._end_node_id}', got '{ns.node_id}'"
                )

    def _build_and_analyze(self) -> None:
        member = Member(
            id=self.id,
            start=self.start.copy(),
            end=self.end.copy(),
            section=self.section,
            material=self.material,
            orientation=self.orientation,
            element_type=self.element_type,
            releases=self.releases,
            constraints=self.constraints,
        )

        # Build frame from member (auto-generates nodes)
        frame = Frame.from_members([member])
        
        # Apply end supports to the auto-generated nodes
        # The frame auto-generates nodes as N0, N1
        frame.nodes["N0"].support = self.support_start
        frame.nodes["N1"].support = self.support_end
        self._start_node_id = "N0"
        self._end_node_id = "N1"
        settings = self.settings if self.settings is not None else FrameAnalysisSettings()
        frame.analyze(self.load_case, settings=settings)
        self._frame = frame

    def analyze(self) -> FrameAnalysisResult:
        return self._frame.analysis_result

    @property
    def analysis_result(self) -> FrameAnalysisResult:
        return self._frame.analysis_result

    @property
    def frame(self) -> Frame:
        """Access the underlying 2-node / 1-member frame model (with cached analysis results)."""
        return self._frame

    def actions(self, points: int = 201):
        """Return a `MemberActionProfile` for this member (continuous across splits)."""
        return self._frame.demand_provider.actions(self.id, points=points)

    def member_demand(self) -> MemberDemand:
        return self._frame.member_demand(self.id)

    def shear(self, axis: Axis, points: int = 100) -> AnalysisResult:
        return self._transverse_analysis(axis=axis, points=points, mode="shear")

    def bending(self, axis: Axis, points: int = 100) -> AnalysisResult:
        return self._transverse_analysis(axis=axis, points=points, mode="bending")

    def axial(self, points: int = 100) -> AnalysisResult:
        return self._axial_torsion_analysis(points=points, mode="axial")

    def torsion(self, points: int = 100) -> AnalysisResult:
        return self._axial_torsion_analysis(points=points, mode="torsion")

    def deflection(self, axis: Axis, points: int = 100) -> Result:
        return self.bending(axis=axis, points=points).displacement

    def von_mises(self, points: int = 201) -> Result:
        """Compute a simple von Mises envelope along the member centerline."""
        profile = self.actions(points=points)

        A = self.section.A
        Iy = self.section.Iy
        Iz = self.section.Iz
        J = self.section.J
        y_max = self.section.y_max
        z_max = self.section.z_max
        r_max = max(abs(y_max), abs(z_max))

        sigma_axial = profile.axial._values / A if A > 0 else np.zeros_like(profile.axial._values)
        sigma_bending_y = np.abs(profile.bending_y._values) * z_max / Iy if Iy > 0 else np.zeros_like(profile.bending_y._values)
        sigma_bending_z = np.abs(profile.bending_z._values) * y_max / Iz if Iz > 0 else np.zeros_like(profile.bending_z._values)
        sigma = np.abs(sigma_axial) + sigma_bending_y + sigma_bending_z

        tau_shear = (np.abs(profile.shear_y._values) + np.abs(profile.shear_z._values)) / A if A > 0 else np.zeros_like(profile.shear_y._values)
        tau_torsion = np.abs(profile.torsion._values) * r_max / J if J > 0 else np.zeros_like(profile.torsion._values)
        tau = tau_shear + tau_torsion

        return Result(profile.axial._x, np.sqrt(sigma**2 + 3.0 * tau**2))

    def check(self, check_module: Any, **kwargs) -> Any:
        """Run a check module on this LoadedMember."""
        if hasattr(check_module, "run"):
            return check_module.run(self, **kwargs)
        if callable(check_module):
            return check_module(self, **kwargs)
        raise ValueError("check_module must be a module with a run() function or a callable.")

    def plot(self, **kwargs) -> None:
        from ..viz.beam_plots import plot_beam_diagram

        plot_beam_diagram(self, **kwargs)

    def plot_results(self, **kwargs) -> None:
        from ..viz.beam_plots import plot_analysis_results

        plot_analysis_results(self, **kwargs)

    def _transverse_analysis(self, axis: Axis, points: int, mode: Literal["shear", "bending"]) -> AnalysisResult:
        profile = self.actions(points=points)

        if axis == "y":
            I = self.section.Iz
            c = self.section.y_max
            action = profile.shear_y if mode == "shear" else profile.bending_z
        else:
            I = self.section.Iy
            c = self.section.z_max
            action = profile.shear_z if mode == "shear" else profile.bending_y

        A = self.section.A
        if mode == "shear":
            stress_values = action._values / A if A > 0 else action._values * 0.0
        else:
            stress_values = action._values * c / I if I > 0 else action._values * 0.0

        disp_values = np.zeros_like(action._values)

        return AnalysisResult(
            Result(action._x, action._values),
            Result(action._x, stress_values),
            Result(action._x, disp_values),
        )

    def _axial_torsion_analysis(self, points: int, mode: Literal["axial", "torsion"]) -> AnalysisResult:
        profile = self.actions(points=points)
        if mode == "axial":
            action = profile.axial
            prop = float(self.section.A)
            stress_factor = 1.0 / prop if prop > 0 else 0.0
        else:
            action = profile.torsion
            prop = float(self.section.J)
            r_max = max(abs(self.section.y_max), abs(self.section.z_max))
            stress_factor = r_max / prop if prop > 0 else 0.0

        stress_values = action._values * stress_factor
        disp_values = np.zeros_like(action._values)

        return AnalysisResult(
            Result(action._x, action._values),
            Result(action._x, stress_values),
            Result(action._x, disp_values),
        )


