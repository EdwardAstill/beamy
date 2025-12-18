# analysis.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Tuple, Optional
import numpy as np

from .frame import Frame
from .loads import FrameLoadCase, MemberPointForce, MemberDistributedForce
from .member import Member


def _build_local_stiffness_matrix(member: Member) -> np.ndarray:
    """
    Build 12×12 local stiffness matrix for a 3D beam element.
    
    DOFs per node (6): [Ux, Uy, Uz, Rx, Ry, Rz]
    Total DOFs: 12 (6 at start node + 6 at end node)
    
    Uses Euler-Bernoulli beam theory:
    - Axial: EA/L
    - Torsion: GJ/L  
    - Bending (y-axis): 12EIz/L³ (shear) + 6EIz/L² (moment coupling)
    - Bending (z-axis): 12EIy/L³ (shear) + 6EIy/L² (moment coupling)
    
    Args:
        member: Member object with section and material properties
        
    Returns:
        12×12 local stiffness matrix
    """
    L = member.length
    E = member.material.E
    G = member.material.G
    A = member.section.A
    Iy = member.section.Iy  # Second moment about y-axis (bending in xz plane)
    Iz = member.section.Iz  # Second moment about z-axis (bending in xy plane)
    J = member.section.J    # Torsion constant
    
    # Initialize matrix
    k = np.zeros((12, 12))
    
    # Axial stiffness (DOFs 0 and 6: Ux at start and end)
    EA_L = E * A / L
    k[0, 0] = EA_L
    k[0, 6] = -EA_L
    k[6, 0] = -EA_L
    k[6, 6] = EA_L
    
    # Torsional stiffness (DOFs 3 and 9: Rx at start and end)
    GJ_L = G * J / L
    k[3, 3] = GJ_L
    k[3, 9] = -GJ_L
    k[9, 3] = -GJ_L
    k[9, 9] = GJ_L
    
    # Bending in xy-plane (about z-axis)
    # DOFs: Uy (1, 7) and Rz (5, 11)
    EIz_L3 = 12 * E * Iz / (L ** 3)
    EIz_L2 = 6 * E * Iz / (L ** 2)
    EIz_L = 4 * E * Iz / L
    EIz_L_half = 2 * E * Iz / L
    
    k[1, 1] = EIz_L3
    k[1, 5] = EIz_L2
    k[1, 7] = -EIz_L3
    k[1, 11] = EIz_L2
    
    k[5, 1] = EIz_L2
    k[5, 5] = EIz_L
    k[5, 7] = -EIz_L2
    k[5, 11] = EIz_L_half
    
    k[7, 1] = -EIz_L3
    k[7, 5] = -EIz_L2
    k[7, 7] = EIz_L3
    k[7, 11] = -EIz_L2
    
    k[11, 1] = EIz_L2
    k[11, 5] = EIz_L_half
    k[11, 7] = -EIz_L2
    k[11, 11] = EIz_L
    
    # Bending in xz-plane (about y-axis)
    # DOFs: Uz (2, 8) and Ry (4, 10)
    EIy_L3 = 12 * E * Iy / (L ** 3)
    EIy_L2 = 6 * E * Iy / (L ** 2)
    EIy_L = 4 * E * Iy / L
    EIy_L_half = 2 * E * Iy / L
    
    k[2, 2] = EIy_L3
    k[2, 4] = -EIy_L2  # Note: sign convention for right-hand rule
    k[2, 8] = -EIy_L3
    k[2, 10] = -EIy_L2
    
    k[4, 2] = -EIy_L2
    k[4, 4] = EIy_L
    k[4, 8] = EIy_L2
    k[4, 10] = EIy_L_half
    
    k[8, 2] = -EIy_L3
    k[8, 4] = EIy_L2
    k[8, 8] = EIy_L3
    k[8, 10] = EIy_L2
    
    k[10, 2] = -EIy_L2
    k[10, 4] = EIy_L_half
    k[10, 8] = EIy_L2
    k[10, 10] = EIy_L
    
    return k


def _build_transformation_matrix_12x12(member: Member) -> np.ndarray:
    """
    Build 12×12 transformation matrix from local to global coordinates.
    
    The 3×3 rotation matrix is expanded to 12×12 by placing it in 4 blocks
    along the diagonal (one for each set of 3 DOFs: translations and rotations
    at start node, then at end node).
    
    Args:
        member: Member object with transformation matrix
        
    Returns:
        12×12 transformation matrix
    """
    T3 = member.transformation_matrix  # 3×3 rotation matrix
    
    # Build 12×12 transformation matrix with 3×3 blocks on diagonal
    T = np.zeros((12, 12))
    
    # Start node translations (DOFs 0-2)
    T[0:3, 0:3] = T3
    # Start node rotations (DOFs 3-5)
    T[3:6, 3:6] = T3
    # End node translations (DOFs 6-8)
    T[6:9, 6:9] = T3
    # End node rotations (DOFs 9-11)
    T[9:12, 9:12] = T3
    
    return T


@dataclass
class LoadedFrame:
    """
    A frame with loads applied, ready for analysis.
    
    On initialization, performs 3D frame analysis using the direct stiffness method.
    
    Attributes:
        frame: Frame geometry (nodes and members)
        loads: Applied loads (nodal and member loads)
    """
    frame: Frame
    loads: FrameLoadCase
    
    # Analysis results (computed in __post_init__)
    _nodal_displacements: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _reactions: Dict[str, np.ndarray] = field(init=False, default_factory=dict)
    _member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(init=False, default_factory=dict)
    
    def __post_init__(self) -> None:
        """Perform frame analysis on initialization."""
        # Validate load references
        self._validate_loads()
        
        # Adjust distributed loads marked as "full length" (end_position = -1)
        self._adjust_full_length_loads()
        
        # Perform 3D frame analysis
        self._analyze_frame()
    
    def _validate_loads(self) -> None:
        """
        Validate that all load references (node IDs and member IDs) exist in the frame.
        
        Raises:
            ValueError: If any load references a non-existent node or member
        """
        # Check nodal forces
        for nf in self.loads.nodal_forces:
            if nf.node_id not in self.frame.nodes:
                raise ValueError(f"NodalForce references unknown node '{nf.node_id}'")
        
        # Check nodal moments
        for nm in self.loads.nodal_moments:
            if nm.node_id not in self.frame.nodes:
                raise ValueError(f"NodalMoment references unknown node '{nm.node_id}'")
        
        # Check member point forces
        for mpf in self.loads.member_point_forces:
            try:
                self.frame.get_member(mpf.member_id)
            except KeyError:
                raise ValueError(f"MemberPointForce references unknown member '{mpf.member_id}'")
        
        # Check member distributed forces
        for mdf in self.loads.member_distributed_forces:
            try:
                self.frame.get_member(mdf.member_id)
            except KeyError:
                raise ValueError(f"MemberDistributedForce references unknown member '{mdf.member_id}'")
    
    def _adjust_full_length_loads(self) -> None:
        """
        Adjust distributed loads marked with end_position=-1 to use actual member length.
        """
        for mdf in self.loads.member_distributed_forces:
            if mdf.end_position == -1.0:
                member = self.frame.get_member(mdf.member_id)
                mdf.end_position = member.length
    
    def _analyze_frame(self) -> None:
        """
        Perform 3D frame analysis using the direct stiffness method.
        
        Steps:
        1. Build global stiffness matrix
        2. Build global load vector (transform member loads to nodal loads)
        3. Apply boundary conditions
        4. Solve for displacements
        5. Compute reactions
        6. Extract member end forces
        """
        # Create node index mapping
        node_ids = sorted(self.frame.nodes.keys())
        node_to_idx = {node_id: i for i, node_id in enumerate(node_ids)}
        n_nodes = len(node_ids)
        n_dofs = 6 * n_nodes
        
        # 1. Build global stiffness matrix
        K_global = np.zeros((n_dofs, n_dofs))
        
        for member in self.frame.members:
            # Build local stiffness matrix
            k_local = _build_local_stiffness_matrix(member)
            
            # Build transformation matrix
            T = _build_transformation_matrix_12x12(member)
            
            # Transform to global coordinates: k_global = T^T @ k_local @ T
            k_global = T.T @ k_local @ T
            
            # Get global DOF indices for this member
            start_idx = node_to_idx[member.start_node_id]
            end_idx = node_to_idx[member.end_node_id]
            
            # DOF indices: [start node 6 DOFs, end node 6 DOFs]
            dof_indices = list(range(start_idx * 6, start_idx * 6 + 6)) + \
                         list(range(end_idx * 6, end_idx * 6 + 6))
            
            # Assemble into global matrix
            for i in range(12):
                for j in range(12):
                    K_global[dof_indices[i], dof_indices[j]] += k_global[i, j]
        
        # 2. Build global load vector
        F_global = self._build_load_vector(node_to_idx, n_dofs)
        
        # 3. Apply boundary conditions
        fixed_dofs = []
        for node_id, node in self.frame.nodes.items():
            if node.support is not None:
                node_idx = node_to_idx[node_id]
                for dof in range(6):
                    if node.support[dof] == '1':
                        fixed_dofs.append(node_idx * 6 + dof)
        
        fixed_dofs = np.array(sorted(fixed_dofs), dtype=int)
        all_dofs = np.arange(n_dofs, dtype=int)
        free_dofs = np.array([d for d in all_dofs if d not in fixed_dofs], dtype=int)
        
        # 4. Solve for displacements
        d_global = np.zeros(n_dofs)
        
        if free_dofs.size > 0:
            K_ff = K_global[np.ix_(free_dofs, free_dofs)]
            F_f = F_global[free_dofs]
            
            # Solve reduced system
            try:
                d_f = np.linalg.solve(K_ff, F_f)
                d_global[free_dofs] = d_f
            except np.linalg.LinAlgError as e:
                raise ValueError(f"Frame is unstable or ill-conditioned: {e}")
        
        # Store nodal displacements
        for node_id in node_ids:
            node_idx = node_to_idx[node_id]
            dof_start = node_idx * 6
            self._nodal_displacements[node_id] = d_global[dof_start:dof_start + 6]
        
        # 5. Compute reactions
        R_global = K_global @ d_global - F_global
        
        for node_id, node in self.frame.nodes.items():
            if node.support is not None:
                node_idx = node_to_idx[node_id]
                dof_start = node_idx * 6
                self._reactions[node_id] = R_global[dof_start:dof_start + 6]
        
        # 6. Extract member end forces (in local coordinates)
        for member in self.frame.members:
            start_idx = node_to_idx[member.start_node_id]
            end_idx = node_to_idx[member.end_node_id]
            
            # Get global displacements at member ends
            d_start_global = d_global[start_idx * 6:(start_idx + 1) * 6]
            d_end_global = d_global[end_idx * 6:(end_idx + 1) * 6]
            d_member_global = np.concatenate([d_start_global, d_end_global])
            
            # Transform to local coordinates
            T = _build_transformation_matrix_12x12(member)
            d_member_local = T @ d_member_global
            
            # Apply local stiffness matrix to get member forces
            k_local = _build_local_stiffness_matrix(member)
            f_member_local = k_local @ d_member_local
            
            # Store as [Nx, Vy, Vz, Tx, My, Mz] at each end
            # Note: DOF ordering is [Ux, Uy, Uz, Rx, Ry, Rz]
            # Force ordering should match: [Fx, Fy, Fz, Mx, My, Mz]
            start_forces = f_member_local[0:6]
            end_forces = f_member_local[6:12]
            
            self._member_end_forces[member.id] = (start_forces, end_forces)
    
    def _build_load_vector(self, node_to_idx: Dict[str, int], n_dofs: int) -> np.ndarray:
        """
        Build global load vector from nodal and member loads.
        
        Member loads are converted to equivalent nodal loads.
        
        Args:
            node_to_idx: Mapping from node ID to node index
            n_dofs: Total number of DOFs
            
        Returns:
            Global load vector
        """
        F = np.zeros(n_dofs)
        
        # Add nodal forces (already in global coordinates)
        for nf in self.loads.nodal_forces:
            node_idx = node_to_idx[nf.node_id]
            F[node_idx * 6:node_idx * 6 + 3] += nf.force
        
        # Add nodal moments (already in global coordinates)
        for nm in self.loads.nodal_moments:
            node_idx = node_to_idx[nm.node_id]
            F[node_idx * 6 + 3:node_idx * 6 + 6] += nm.moment
        
        # Convert member point forces to equivalent nodal loads
        for mpf in self.loads.member_point_forces:
            self._add_member_point_force_to_vector(mpf, node_to_idx, F)
        
        # Convert member distributed forces to equivalent nodal loads
        for mdf in self.loads.member_distributed_forces:
            self._add_member_distributed_force_to_vector(mdf, node_to_idx, F)
        
        return F
    
    def _add_member_point_force_to_vector(
        self,
        mpf: MemberPointForce,
        node_to_idx: Dict[str, int],
        F: np.ndarray
    ) -> None:
        """
        Add a member point force to the global load vector as equivalent nodal loads.
        
        Uses linear interpolation (shape functions) to distribute the force to nodes.
        """
        member = self.frame.get_member(mpf.member_id)
        
        # Get position along member
        if mpf.position_type == "relative":
            x = mpf.position * member.length
        else:
            x = mpf.position
        
        # Transform force to global coordinates if needed
        if mpf.coords == "local":
            T3 = member.transformation_matrix
            force_global = T3.T @ mpf.force
        else:
            force_global = mpf.force
        
        # Linear shape functions: N1 = 1 - x/L, N2 = x/L
        L = member.length
        N1 = 1.0 - x / L
        N2 = x / L
        
        # Distribute force to nodes
        start_idx = node_to_idx[member.start_node_id]
        end_idx = node_to_idx[member.end_node_id]
        
        F[start_idx * 6:start_idx * 6 + 3] += N1 * force_global
        F[end_idx * 6:end_idx * 6 + 3] += N2 * force_global
    
    def _add_member_distributed_force_to_vector(
        self,
        mdf: MemberDistributedForce,
        node_to_idx: Dict[str, int],
        F: np.ndarray
    ) -> None:
        """
        Add a member distributed force to the global load vector.
        
        Uses equivalent nodal loads based on distributed load integration.
        For simplicity, discretizes the distributed load into point loads.
        """
        member = self.frame.get_member(mdf.member_id)
        
        # Discretize distributed load into equivalent point loads
        n_points = 11  # Number of integration points
        x_positions = np.linspace(mdf.start_position, mdf.end_position, n_points)
        dx = (mdf.end_position - mdf.start_position) / (n_points - 1) if n_points > 1 else 0
        
        for i, x in enumerate(x_positions):
            # Interpolate force at this position
            t = (x - mdf.start_position) / (mdf.end_position - mdf.start_position) if mdf.end_position > mdf.start_position else 0
            force_per_length = (1 - t) * mdf.start_force + t * mdf.end_force
            
            # Convert to point force
            point_force = force_per_length * dx
            
            # Transform to global if needed
            if mdf.coords == "local":
                T3 = member.transformation_matrix
                point_force = T3.T @ point_force
            
            # Distribute to nodes using shape functions
            L = member.length
            N1 = 1.0 - x / L
            N2 = x / L
            
            start_idx = node_to_idx[member.start_node_id]
            end_idx = node_to_idx[member.end_node_id]
            
            F[start_idx * 6:start_idx * 6 + 3] += N1 * point_force
            F[end_idx * 6:end_idx * 6 + 3] += N2 * point_force
    
    @property
    def nodal_displacements(self) -> Dict[str, np.ndarray]:
        """
        Dictionary mapping node ID to displacement vector [UX, UY, UZ, RX, RY, RZ].
        Displacements are in global coordinates.
        """
        return self._nodal_displacements
    
    @property
    def reactions(self) -> Dict[str, np.ndarray]:
        """
        Dictionary mapping supported node ID to reaction vector [FX, FY, FZ, MX, MY, MZ].
        Reactions are in global coordinates.
        """
        return self._reactions
    
    @property
    def member_end_forces(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Dictionary mapping member ID to (start_forces, end_forces).
        Forces are in local member coordinates [Fx, Fy, Fz, Mx, My, Mz].
        """
        return self._member_end_forces
    
    def to_loaded_beams(self) -> Dict[str, 'LoadedBeam']:
        """
        Convert frame analysis to individual LoadedBeam objects.
        
        Each member is converted to a LoadedBeam with:
        - Member end forces applied as support reactions
        - Any intermediate member loads from the original FrameLoadCase
        - Local coordinate system orientation
        
        Returns:
            Dictionary mapping member ID to LoadedBeam
        """
        from ..setup.beam import Beam1D, Support
        from ..setup.loads import LoadCase, PointForce, DistributedForce
        from ..analysis.analysis import LoadedBeam
        
        loaded_beams = {}
        
        for member in self.frame.members:
            # Get member end forces (in local coordinates)
            start_forces, end_forces = self._member_end_forces[member.id]
            
            # Create a Beam1D for this member
            # Supports at both ends with reactions
            beam = Beam1D(
                L=member.length,
                material=member.material,
                section=member.section,
                supports=[
                    Support(x=0.0, type="111111"),  # Fixed at start
                    Support(x=member.length, type="111111")  # Fixed at end
                ]
            )
            
            # Create a LoadCase with member end forces applied as external loads
            # (opposite sign to balance the reactions)
            load_case = LoadCase(name=f"Member {member.id} loads")
            
            # Add start end forces as point loads (negated to represent applied loads)
            # Format: [Fx, Fy, Fz, Mx, My, Mz]
            load_case.add_point_force(PointForce(
                point=np.array([0.0, 0.0, 0.0]),
                force=-start_forces[0:3]  # [Fx, Fy, Fz] at start
            ))
            
            load_case.add_point_force(PointForce(
                point=np.array([member.length, 0.0, 0.0]),
                force=-end_forces[0:3]  # [Fx, Fy, Fz] at end
            ))
            
            # Add any intermediate member loads from the original load case
            for mpf in self.loads.member_point_forces:
                if mpf.member_id == member.id:
                    # Get position
                    if mpf.position_type == "relative":
                        x = mpf.position * member.length
                    else:
                        x = mpf.position
                    
                    # Transform force to local if needed
                    if mpf.coords == "global":
                        T3 = member.transformation_matrix
                        force_local = T3 @ mpf.force
                    else:
                        force_local = mpf.force
                    
                    load_case.add_point_force(PointForce(
                        point=np.array([x, 0.0, 0.0]),
                        force=force_local
                    ))
            
            # Add distributed loads
            for mdf in self.loads.member_distributed_forces:
                if mdf.member_id == member.id:
                    # Transform forces to local if needed
                    if mdf.coords == "global":
                        T3 = member.transformation_matrix
                        start_force_local = T3 @ mdf.start_force
                        end_force_local = T3 @ mdf.end_force
                    else:
                        start_force_local = mdf.start_force
                        end_force_local = mdf.end_force
                    
                    load_case.add_distributed_force(DistributedForce(
                        start_position=np.array([mdf.start_position, 0.0, 0.0]),
                        end_position=np.array([mdf.end_position, 0.0, 0.0]),
                        start_force=start_force_local,
                        end_force=end_force_local
                    ))
            
            # Create LoadedBeam
            loaded_beam = LoadedBeam(beam, load_case)
            loaded_beams[member.id] = loaded_beam
        
        return loaded_beams
    
    def get_member_results(self, member_id: str) -> 'MemberResults':
        """
        Get detailed analysis results for a specific member.
        
        Args:
            member_id: Member identifier
            
        Returns:
            MemberResults object with detailed force/stress/displacement distributions
        """
        member = self.frame.get_member(member_id)
        loaded_beams = self.to_loaded_beams()
        loaded_beam = loaded_beams[member_id]
        
        return MemberResults(member=member, loaded_beam=loaded_beam)
    
    def plot(
        self,
        show_loads: bool = True,
        show_reactions: bool = True,
        show_member_ids: bool = True,
        show_node_ids: bool = True,
        deformed: bool = False,
        scale_factor: float = 1.0,
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot the frame in 3D wireframe style.
        
        Args:
            show_loads: Display applied load arrows
            show_reactions: Display reaction arrows at supports
            show_member_ids: Label members with their IDs
            show_node_ids: Label nodes with their IDs
            deformed: If True, show deformed shape overlay
            scale_factor: Scale factor for deformed shape visualization
            save_path: Path to save the plot (supports .svg)
        """
        from .plotter import plot_frame
        plot_frame(self, show_loads, show_reactions, show_member_ids, show_node_ids,
                  deformed, scale_factor, save_path)
    
    def plot_deflection(
        self,
        scale_factor: float = 1.0,
        points_per_member: int = 20,
        colormap: str = "viridis",
        show_undeformed: bool = True,
        show_colorbar: bool = True,
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot the deformed frame shape in 3D wireframe, colored by displacement magnitude.
        
        Args:
            scale_factor: Multiplier for displacement visualization
            points_per_member: Number of interpolation points along each member
            colormap: Matplotlib colormap name
            show_undeformed: Show original geometry as faint dashed lines
            show_colorbar: Display colorbar with displacement units
            save_path: Path to save the plot (supports .svg)
        """
        from .plotter import plot_deflection
        plot_deflection(self, scale_factor, points_per_member, colormap,
                       show_undeformed, show_colorbar, save_path)
    
    def plot_von_mises(
        self,
        points_per_member: int = 20,
        colormap: str = "jet",
        show_colorbar: bool = True,
        stress_limits: Optional[Tuple[float, float]] = None,
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot the frame in 3D wireframe, colored by Von Mises stress.
        
        Args:
            points_per_member: Number of interpolation points along each member
            colormap: Matplotlib colormap name
            show_colorbar: Display colorbar with stress units
            stress_limits: Optional (min, max) to fix colorbar range
            save_path: Path to save the plot (supports .svg)
        """
        from .plotter import plot_von_mises
        plot_von_mises(self, points_per_member, colormap, show_colorbar,
                      stress_limits, save_path)
    
    def plot_results(
        self,
        result_type: str = "von_mises",
        deformed: bool = True,
        scale_factor: float = 1.0,
        points_per_member: int = 20,
        colormap: str = "jet",
        show_undeformed: bool = True,
        show_colorbar: bool = True,
        show_node_ids: bool = False,
        show_member_ids: bool = False,
        value_limits: Optional[Tuple[float, float]] = None,
        save_path: Optional[str] = None
    ) -> None:
        """
        Unified 3D wireframe plot showing deformed shape colored by analysis results.
        
        Args:
            result_type: Type of result - "von_mises" or "deflection"
            deformed: Plot on deformed geometry
            scale_factor: Displacement scale factor
            points_per_member: Number of interpolation points per member
            colormap: Matplotlib colormap name
            show_undeformed: Show original geometry as reference
            show_colorbar: Display colorbar legend
            show_node_ids: Label nodes
            show_member_ids: Label members
            value_limits: Optional (min, max) to fix color range
            save_path: Path to save the plot
        """
        from .plotter import plot_results
        plot_results(self, result_type, deformed, scale_factor, points_per_member,
                    colormap, show_undeformed, show_colorbar, show_node_ids,
                    show_member_ids, value_limits, save_path)
    
    def plot_member_diagrams(
        self,
        member_id: str,
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot internal force diagrams (N, V, M, T) for a specific member.
        
        Args:
            member_id: Member identifier
            save_path: Path to save the plot
        """
        from .plotter import plot_member_diagrams
        plot_member_diagrams(self, member_id, save_path)


@dataclass
class MemberResults:
    """
    Detailed analysis results for an individual member.
    
    Provides access to force, stress, and displacement distributions along the member
    by wrapping a LoadedBeam object.
    
    Attributes:
        member: Member object
        loaded_beam: LoadedBeam object with analysis results (private)
    """
    member: Member
    loaded_beam: 'LoadedBeam'
    
    @property
    def axial(self):
        """Axial force distribution: N(x), sigma_axial(x), u(x)."""
        return self.loaded_beam.axial()
    
    @property
    def shear_y(self):
        """Shear force in y-direction: Vy(x), tau_y(x), v(x)."""
        return self.loaded_beam.shear("y")
    
    @property
    def shear_z(self):
        """Shear force in z-direction: Vz(x), tau_z(x), w(x)."""
        return self.loaded_beam.shear("z")
    
    @property
    def bending_y(self):
        """Bending about y-axis: My(x), sigma_y(x), theta_y(x)."""
        return self.loaded_beam.bending("y")
    
    @property
    def bending_z(self):
        """Bending about z-axis: Mz(x), sigma_z(x), theta_z(x)."""
        return self.loaded_beam.bending("z")
    
    @property
    def torsion(self):
        """Torsional moment: T(x), tau_torsion(x), phi(x)."""
        return self.loaded_beam.torsion()
    
    @property
    def von_mises(self):
        """Combined Von Mises stress distribution."""
        return self.loaded_beam.von_mises()
