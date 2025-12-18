from __future__ import annotations
from typing import Dict, List, Tuple
import numpy as np

from ..core.math import build_local_stiffness_matrix, build_transformation_matrix_12x12
from .frame import Frame
from .member import Member

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

def analyze_frame_geometry(frame: Frame, node_to_idx: Dict[str, int]) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Build global stiffness matrix for a frame.
    """
    n_nodes = len(frame.nodes)
    n_dofs = 6 * n_nodes
    K_global = np.zeros((n_dofs, n_dofs))
    member_matrices = {}

    for member in frame.members:
        # Build matrices
        k_local = build_local_stiffness_matrix(
            member.length, member.material.E, member.material.G,
            member.section.A, member.section.Iy, member.section.Iz, member.section.J
        )
        
        # Apply member end releases if specified
        if member.releases:
            k_local = apply_releases(k_local, member.releases)
        
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

