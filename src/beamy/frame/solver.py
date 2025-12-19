from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import numpy as np

from ..core.math import build_local_stiffness_matrix, build_transformation_matrix_12x12
from .frame import Frame
from .member import Member


MemberStiffnessMatrices = Dict[str, Tuple[np.ndarray, np.ndarray]]


@dataclass(frozen=True)
class ElementStiffnessScales:
    """Per-element stiffness scaling.

    These scale the section/material properties used to build the element stiffness.
    Defaults are 1.0 (no modification).
    """

    E: float = 1.0
    G: float = 1.0
    A: float = 1.0
    Iy: float = 1.0
    Iz: float = 1.0
    J: float = 1.0

def build_local_truss_stiffness_matrix(L: float, E: float, A: float, axial_scale: float = 1.0) -> np.ndarray:
    """
    Build a 12x12 local stiffness matrix for an axial-only element (truss/cable) embedded
    in the 3D frame DOF set [UX UY UZ RX RY RZ] at each node.

    Only the local axial translation DOFs (start UX = 0, end UX = 6) carry stiffness.
    All other DOFs have zero stiffness (they are not coupled by this element).

    axial_scale can be used to "slacken" a cable (e.g. in compression).
    """
    k = np.zeros((12, 12))
    if L <= 0.0:
        raise ValueError("L must be positive")
    EA_L = (E * A / L) * axial_scale
    k[0, 0] = EA_L
    k[0, 6] = -EA_L
    k[6, 0] = -EA_L
    k[6, 6] = EA_L
    return k


def apply_releases(k_local: np.ndarray, releases: str) -> np.ndarray:
    """
    Apply member end releases to local stiffness matrix.
    
    Releases is a 12-character string of '0' and '1':
    - First 6 chars: start node [UX, UY, UZ, RX, RY, RZ]
    - Last 6 chars: end node [UX, UY, UZ, RX, RY, RZ]
    - '0' = connected (rigid), '1' = released (free)
    
    When a DOF is released, that DOF cannot transfer forces/moments.
    We zero out the corresponding rows and columns in the stiffness matrix.
    """
    k_modified = k_local.copy()
    
    for i, release_char in enumerate(releases):
        if release_char == '1':
            # Zero out row and column for released DOF
            k_modified[i, :] = 0.0
            k_modified[:, i] = 0.0
            # Add small diagonal value to prevent singularity
            k_modified[i, i] = 1e-6
    
    return k_modified

def analyze_frame_geometry(
    frame: Frame,
    node_to_idx: Dict[str, int],
    member_axial_scales: Optional[Dict[str, float]] = None,
    member_stiffness_scales: Optional[Dict[str, ElementStiffnessScales]] = None,
) -> Tuple[np.ndarray, MemberStiffnessMatrices]:
    """Assemble the global stiffness matrix for a frame.

    Returns:
        (K_global, member_matrices)

    where member_matrices[member_id] = (k_local, T_12x12).
    """
    n_nodes = len(frame.nodes)
    n_dofs = 6 * n_nodes
    K_global = np.zeros((n_dofs, n_dofs))
    member_matrices: MemberStiffnessMatrices = {}

    for member in frame.members:
        axial_scale = 1.0
        if member_axial_scales is not None and member.id in member_axial_scales:
            axial_scale = float(member_axial_scales[member.id])

        stiff = ElementStiffnessScales()
        if member_stiffness_scales is not None and member.id in member_stiffness_scales:
            stiff = member_stiffness_scales[member.id]

        # Build matrices (beam vs truss/cable)
        if member.element_type == "beam":
            k_local = build_local_stiffness_matrix(
                member.length,
                member.material.E * stiff.E,
                member.material.G * stiff.G,
                member.section.A * stiff.A,
                member.section.Iy * stiff.Iy,
                member.section.Iz * stiff.Iz,
                member.section.J * stiff.J,
            )
            # Apply member end releases only for beam elements
            if member.releases:
                k_local = apply_releases(k_local, member.releases)
        else:
            # truss/cable: axial-only
            k_local = build_local_truss_stiffness_matrix(
                L=member.length,
                E=member.material.E,
                A=member.section.A,
                axial_scale=axial_scale,
            )
        
        T = build_transformation_matrix_12x12(member.transformation_matrix)
        k_global = T.T @ k_local @ T
        
        # Assemble
        s_idx, e_idx = node_to_idx[member.start_node_id], node_to_idx[member.end_node_id]
        dof_indices = list(range(s_idx * 6, s_idx * 6 + 6)) + list(range(e_idx * 6, e_idx * 6 + 6))
        for i in range(12):
            for j in range(12):
                K_global[dof_indices[i], dof_indices[j]] += k_global[i, j]
        
        member_matrices[member.id] = (k_local, T)
        
    return K_global, member_matrices


def assemble_global_stiffness(
    frame: Frame,
    node_to_idx: Dict[str, int],
    member_axial_scales: Optional[Dict[str, float]] = None,
    member_stiffness_scales: Optional[Dict[str, ElementStiffnessScales]] = None,
) -> Tuple[np.ndarray, MemberStiffnessMatrices]:
    """Public wrapper with a clearer name.

    This exists so the higher-level analysis code can follow a stable
    assemble → solve → recover structure.
    """

    return analyze_frame_geometry(
        frame,
        node_to_idx,
        member_axial_scales=member_axial_scales,
        member_stiffness_scales=member_stiffness_scales,
    )

def solve_displacements(K_global: np.ndarray, F_global: np.ndarray, fixed_dofs: List[int]) -> np.ndarray:
    """
    Solve for global displacements.
    """
    n_dofs = len(F_global)
    d_global = np.zeros(n_dofs)
    free_dofs = [i for i in range(n_dofs) if i not in fixed_dofs]
    
    if free_dofs:
        K_ff = K_global[np.ix_(free_dofs, free_dofs)]
        F_f = F_global[free_dofs]
        try:
            d_global[free_dofs] = np.linalg.solve(K_ff, F_f)
        except np.linalg.LinAlgError:
            raise ValueError("Frame is unstable or ill-conditioned.")
            
    return d_global


def recover_member_end_forces(
    frame: Frame,
    node_to_idx: Dict[str, int],
    d_global: np.ndarray,
    member_matrices: MemberStiffnessMatrices,
) -> Tuple[Dict[str, np.ndarray], Dict[str, Tuple[np.ndarray, np.ndarray]]]:
    """Recover nodal displacements and member end forces from a solved displacement vector.

    Returns:
        (nodal_displacements, member_end_forces)

    Nodal displacements are per-node 6-vectors in global axes.
    Member end forces are stored in *local* axes as two 6-vectors:
        start = [Fx, Fy, Fz, Mx, My, Mz]
        end   = [Fx, Fy, Fz, Mx, My, Mz]
    """

    node_ids = sorted(frame.nodes.keys())
    nodal_displacements: Dict[str, np.ndarray] = {
        nid: d_global[node_to_idx[nid] * 6 : node_to_idx[nid] * 6 + 6] for nid in node_ids
    }

    member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for member in frame.members:
        d_m_g = np.concatenate(
            [nodal_displacements[member.start_node_id], nodal_displacements[member.end_node_id]]
        )
        k_l, T = member_matrices[member.id]
        f_m_l = k_l @ (T @ d_m_g)
        member_end_forces[member.id] = (f_m_l[0:6], f_m_l[6:12])

    return nodal_displacements, member_end_forces


def build_local_geometric_stiffness_matrix(L: float, N: float) -> np.ndarray:
    """Build a 12x12 local geometric stiffness matrix for a 3D beam element.

    Conventions used by beamy:
        - DOF order per node: [Ux, Uy, Uz, Rx, Ry, Rz]
        - Internal axial force sign: N > 0 is tension, N < 0 is compression.

    Behavior:
        - Compression (N < 0) is destabilizing (negative contribution to lateral stiffness).
        - Tension (N > 0) is stabilizing.

    Notes:
        - This implements the classic beam-column geometric stiffness for bending
          in the (x,y) and (x,z) planes.
        - Torsional geometric stiffness is not included.
    """

    if L <= 0.0:
        raise ValueError("L must be positive")

    kg = np.zeros((12, 12))
    c = float(N) / (30.0 * float(L))

    # Standard 2D geometric stiffness for DOFs [v1, theta1, v2, theta2].
    L2 = float(L) ** 2
    k2 = np.array(
        [
            [36.0, 3.0 * L, -36.0, 3.0 * L],
            [3.0 * L, 4.0 * L2, -3.0 * L, -1.0 * L2],
            [-36.0, -3.0 * L, 36.0, -3.0 * L],
            [3.0 * L, -1.0 * L2, -3.0 * L, 4.0 * L2],
        ],
        dtype=float,
    )
    k2 *= c

    # Bending in xy-plane (about local z): Uy + Rz match the standard sign convention.
    dofs_y = [1, 5, 7, 11]
    for i, gi in enumerate(dofs_y):
        for j, gj in enumerate(dofs_y):
            kg[gi, gj] += k2[i, j]

    # Bending in xz-plane (about local y): Uz + Ry in beamy has opposite sign vs standard theta.
    # Map theta_standard = -Ry_beamy using similarity transform S.
    S = np.diag([1.0, -1.0, 1.0, -1.0])
    k2_ry = S @ k2 @ S
    dofs_z = [2, 4, 8, 10]
    for i, gi in enumerate(dofs_z):
        for j, gj in enumerate(dofs_z):
            kg[gi, gj] += k2_ry[i, j]

    return kg


def assemble_geometric_stiffness(
    frame: Frame,
    node_to_idx: Dict[str, int],
    member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]],
    member_matrices: Optional[MemberStiffnessMatrices] = None,
) -> np.ndarray:
    """Assemble global geometric stiffness matrix Kg from current member axial forces.

    Uses axial force recovered from member end forces with beamy's convention:
        N = -start_Fx  (tension positive)
    """

    n_nodes = len(frame.nodes)
    n_dofs = 6 * n_nodes
    Kg_global = np.zeros((n_dofs, n_dofs))

    for member in frame.members:
        if member.element_type != "beam":
            continue

        if member.id not in member_end_forces:
            continue
        start_f, _end_f = member_end_forces[member.id]
        N = -float(start_f[0])  # tension positive

        kg_local = build_local_geometric_stiffness_matrix(member.length, N)

        if member.releases:
            kg_local = apply_releases(kg_local, member.releases)

        if member_matrices is not None and member.id in member_matrices:
            _k_local, T = member_matrices[member.id]
        else:
            T = build_transformation_matrix_12x12(member.transformation_matrix)

        kg_global = T.T @ kg_local @ T

        s_idx, e_idx = node_to_idx[member.start_node_id], node_to_idx[member.end_node_id]
        dof_indices = list(range(s_idx * 6, s_idx * 6 + 6)) + list(range(e_idx * 6, e_idx * 6 + 6))
        for i in range(12):
            for j in range(12):
                Kg_global[dof_indices[i], dof_indices[j]] += kg_global[i, j]

    return Kg_global

