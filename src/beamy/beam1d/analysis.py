"""
1D Beam Analysis using Frame Backend (Strategy A).

LoadedMember now delegates to a single-member frame internally for consistent
results with full frame analysis. This eliminates divergence between 1D and 3D
analysis paths.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Any
import numpy as np

from .beam import Beam1D
from ..core.results import Result, AnalysisResult
from ..core.loads import LoadCase


@dataclass
class LoadedMember:
    """1D Beam analysis using frame backend (Strategy A).
    
    This builds a single-member frame internally and delegates to the frame solver,
    ensuring consistent results between 1D and 3D analysis paths.
    
    Args:
        beam: Beam1D object defining geometry, material, section, and supports
        loads: LoadCase containing point forces, moments, and distributed loads
        settings: Optional FrameAnalysisSettings for advanced control
    """
    beam: Beam1D
    loads: LoadCase
    settings: Optional[Any] = None
    
    # Frame analysis state
    _frame_analysis: Any = field(init=False, repr=False)
    _member_id: str = field(init=False, default="M1", repr=False)

    def __post_init__(self):
        """Build single-member frame and analyze."""
        self._initialize_frame_backend()
    
    def _initialize_frame_backend(self):
        """Initialize using frame analysis delegation (Strategy A)."""
        # Import here to avoid circular dependency
        from ..frame.frame import Frame
        from ..frame.node import Node
        from ..frame.member import Member
        from ..core.loads import FrameLoadCase
        from ..frame.analysis import LoadedFrame, FrameAnalysisSettings
        
        # Build single-member frame
        # Collect all x-positions where we need nodes: supports + load locations
        node_positions = {0.0, self.beam.L}
        
        # Support positions
        for support in self.beam.supports:
            node_positions.add(support.x)
        
        # Point force positions
        for pf in self.loads.point_forces:
            node_positions.add(float(pf.point[0]))
        
        # Moment positions
        for mom in self.loads.moments:
            node_positions.add(mom.x)
        
        # Distributed force endpoints
        for df in self.loads.dist_forces:
            node_positions.add(float(df.start_position[0]))
            node_positions.add(float(df.end_position[0]))
        
        # Sort positions
        node_positions_sorted = sorted(node_positions)
        
        # Create nodes with support conditions
        nodes = []
        support_by_x = {s.x: s for s in self.beam.supports}
        
        for i, x in enumerate(node_positions_sorted):
            node_id = f"N{i}"
            position = np.array([x, 0.0, 0.0])
            
            # Apply support if this position has one
            support_str = None
            if x in support_by_x:
                support_str = support_by_x[x].type
            
            nodes.append(Node(id=node_id, position=position, support=support_str))
        
        # Create single member connecting first and last node
        # For now, use the full member - later we might split for intermediate supports
        member = Member(
            id=self._member_id,
            start_node_id=nodes[0].id,
            end_node_id=nodes[-1].id,
            section=self.beam.section,
            material=self.beam.material,
            orientation=np.array([0.0, 0.0, 1.0]),  # Default orientation
            element_type="beam",
        )
        
        # Build frame
        frame = Frame.from_nodes_and_members(nodes, [member])
        
        # Convert loads
        frame_loads = self._convert_loads_to_frame(frame)
        
        # Analyze
        analysis_settings = self.settings if self.settings else FrameAnalysisSettings()
        self._frame_analysis = LoadedFrame(frame=frame, loads=frame_loads, settings=analysis_settings)
        
        # For backward compatibility, populate all_loads from frame results
        self.all_loads = []  # Could reconstruct if needed
    
    def _convert_loads_to_frame(self, frame) -> Any:
        """Convert 1D LoadCase to FrameLoadCase for the single-member frame."""
        from ..core.loads import FrameLoadCase
        
        frame_loads = FrameLoadCase(name=self.loads.name)
        
        # Convert point forces to nodal forces
        # For frame-based backend, point forces need to be applied at nodes
        # We'll use member point forces instead to allow intermediate loading
        for pf in self.loads.point_forces:
            x = float(pf.point[0])
            force_local = pf.force  # [Fx, Fy, Fz]
            frame_loads.add_member_point_force(
                member_id=self._member_id,
                position=x,
                force=force_local,
                coords="local",
                position_type="absolute"
            )
        
        # Convert moments to member point moments
        for mom in self.loads.moments:
            moment_local = mom.moment  # [Mx, My, Mz]
            frame_loads.add_member_point_moment(
                member_id=self._member_id,
                position=mom.x,
                moment=moment_local,
                coords="local",
                position_type="absolute"
            )
        
        # Convert distributed forces
        for df in self.loads.dist_forces:
            start_x = float(df.start_position[0])
            end_x = float(df.end_position[0])
            start_force = df.start_force  # [wx, wy, wz]
            end_force = df.end_force
            
            frame_loads.add_member_distributed_force(
                member_id=self._member_id,
                start_position=start_x,
                end_position=end_x,
                start_force=start_force,
                end_force=end_force,
                coords="local"
            )
        
        return frame_loads

    def shear(self, axis, points=100):
        """Get shear force results along beam."""
        return self._transverse_analysis(axis, points, "shear")
    
    def bending(self, axis, points=100):
        """Get bending moment results along beam."""
        return self._transverse_analysis(axis, points, "bending")
    
    def axial(self, points=100):
        """Get axial force results along beam."""
        return self._axial_torsion_analysis(points, "axial")
    
    def torsion(self, points=100):
        """Get torsion results along beam."""
        return self._axial_torsion_analysis(points, "torsion")
    
    def deflection(self, axis, points=100):
        """Get deflection along beam (from bending analysis)."""
        return self.bending(axis, points).displacement

    def von_mises(self, points=100):
        """Compute von Mises stress along beam."""
        profile = self._frame_analysis.demand_provider.actions(self._member_id, points=points)
        
        # Section properties
        A = self.beam.section.A
        Iy = self.beam.section.Iy
        Iz = self.beam.section.Iz
        J = self.beam.section.J
        y_max = self.beam.section.y_max
        z_max = self.beam.section.z_max
        r_max = max(abs(y_max), abs(z_max))
        
        # Compute stresses
        sigma_axial = profile.axial._values / A if A > 0 else np.zeros_like(profile.axial._values)
        sigma_bending_y = np.abs(profile.bending_y._values) * z_max / Iy if Iy > 0 else np.zeros_like(profile.bending_y._values)
        sigma_bending_z = np.abs(profile.bending_z._values) * y_max / Iz if Iz > 0 else np.zeros_like(profile.bending_z._values)
        sigma = np.abs(sigma_axial) + sigma_bending_y + sigma_bending_z
        
        tau_shear = (np.abs(profile.shear_y._values) + np.abs(profile.shear_z._values)) / A if A > 0 else np.zeros_like(profile.shear_y._values)
        tau_torsion = np.abs(profile.torsion._values) * r_max / J if J > 0 else np.zeros_like(profile.torsion._values)
        tau = tau_shear + tau_torsion
        
        return Result(profile.axial._x, np.sqrt(sigma**2 + 3 * tau**2))

    def check(self, check_module: Any, **kwargs) -> Any:
        """
        Run a check module on this loaded beam.
        The check module must have a run(loaded_beam, **kwargs) function.
        """
        if hasattr(check_module, "run"):
            return check_module.run(self, **kwargs)
        if callable(check_module):
            return check_module(self, **kwargs)
        raise ValueError("check_module must be a module with a run() function or a callable.")

    def check_aisc_chapter_f(self, length_unit, force_unit) -> Any:
        from ..checks import aisc_9
        return self.check(aisc_9, length_unit=length_unit, force_unit=force_unit)

    def check_aisc_chapter_e(
        self,
        length_unit,
        force_unit,
        *,
        K: float = 1.0,
        unbraced_positions: Any = None,
        compression_sign: str = "negative",
    ) -> Any:
        from ..checks import aisc_9
        return self.check(
            aisc_9,
            chapter="e",
            length_unit=length_unit,
            force_unit=force_unit,
            K=K,
            unbraced_positions=unbraced_positions,
            compression_sign=compression_sign,
        )

    def check_aisc_chapter_h(
        self,
        length_unit,
        force_unit,
        *,
        Ky: float = 1.0,
        Kz: float = 1.0,
        unbraced_positions_y: Any = None,
        unbraced_positions_z: Any = None,
        Cmx: Any = None,
        Cmy: Any = None,
        frame_type: str = "braced",
        has_transverse_loading: bool = True,
        compression_sign: str = "negative",
    ) -> Any:
        from ..checks import aisc_9
        return self.check(
            aisc_9,
            chapter="h",
            length_unit=length_unit,
            force_unit=force_unit,
            Ky=Ky,
            Kz=Kz,
            unbraced_positions_y=unbraced_positions_y,
            unbraced_positions_z=unbraced_positions_z,
            Cmx=Cmx,
            Cmy=Cmy,
            frame_type=frame_type,
            has_transverse_loading=has_transverse_loading,
            compression_sign=compression_sign,
        )

    def plot(self, **kwargs):
        from ..viz.beam_plots import plot_beam_diagram
        plot_beam_diagram(self, **kwargs)

    def plot_results(self, **kwargs):
        from ..viz.beam_plots import plot_analysis_results
        plot_analysis_results(self, **kwargs)

    def _transverse_analysis(self, axis, points, mode):
        """Get transverse (shear/bending) results from frame analysis."""
        profile = self._frame_analysis.demand_provider.actions(self._member_id, points=points)
        
        if axis == "y":
            I = self.beam.section.Iz
            c = self.beam.section.y_max
            action = profile.shear_y if mode == "shear" else profile.bending_z
        else:  # axis == "z"
            I = self.beam.section.Iy
            c = self.beam.section.z_max
            action = profile.shear_z if mode == "shear" else profile.bending_y
        
        A = self.beam.section.A
        if mode == "shear":
            stress_values = action._values / A if A > 0 else action._values * 0
        else:  # bending
            stress_values = action._values * c / I if I > 0 else action._values * 0
        
        # Displacement not yet fully supported from frame backend
        disp_values = np.zeros_like(action._values)
        
        return AnalysisResult(
            Result(action._x, action._values),
            Result(action._x, stress_values),
            Result(action._x, disp_values)
        )
    
    def _axial_torsion_analysis(self, points, mode):
        """Get axial/torsion results from frame analysis."""
        profile = self._frame_analysis.demand_provider.actions(self._member_id, points=points)
        
        if mode == "axial":
            action = profile.axial
            prop = self.beam.section.A
            stress_factor = 1.0 / prop if prop > 0 else 0.0
        else:  # torsion
            action = profile.torsion
            prop = self.beam.section.J
            r_max = max(abs(self.beam.section.y_max), abs(self.beam.section.z_max))
            stress_factor = r_max / prop if prop > 0 else 0.0
        
        stress_values = action._values * stress_factor
        disp_values = np.zeros_like(action._values)
        
        return AnalysisResult(
            Result(action._x, action._values),
            Result(action._x, stress_values),
            Result(action._x, disp_values)
        )
