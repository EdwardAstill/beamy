import numpy as np

def build_transformation_matrix_12x12(T3: np.ndarray) -> np.ndarray:
    """
    Build 12x12 transformation matrix from 3x3 rotation matrix.
    
    Args:
        T3: 3x3 rotation matrix
        
    Returns:
        12x12 transformation matrix
    """
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

def build_local_stiffness_matrix(
    L: float, E: float, G: float, A: float, Iy: float, Iz: float, J: float
) -> np.ndarray:
    """
    Build 12x12 local stiffness matrix for a 3D beam element.
    """
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
    EIy_L3 = 12 * E * Iy / (L ** 3)
    EIy_L2 = 6 * E * Iy / (L ** 2)
    EIy_L = 4 * E * Iy / L
    EIy_L_half = 2 * E * Iy / L
    
    k[2, 2] = EIy_L3
    k[2, 4] = -EIy_L2
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

