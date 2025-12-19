# Finite Element Method (FEM) Implementation in Beamy

## Overview
Beamy uses the **Displacement Method** (also called stiffness method) to solve structural problems. This is the standard modern FEM approach where displacements are the primary unknowns.

---

## Matrix Assembly & Solution Method

### System Solver
- **Direct Linear Solver**: `np.linalg.solve()` using **LU Decomposition**
  - Applies to the reduced system (after enforcing boundary conditions)
  - Equations: `K_ff × d_f = F_f` where subscript `f` denotes free DOFs
  - Not an iterative method like conjugate gradient

### Assembly Process
1. **Direct Assembly** (Direct Stiffness Method)
   - For each element: compute local stiffness matrix `k_local` → transform to global coordinates → add contributions to global `K_global`
   - Assembly loop: `K_global[dof_indices[i], dof_indices[j]] += k_global[i, j]`
   - Global stiffness matrix is symmetric and positive semi-definite

2. **DOF Numbering Convention** (6 DOF per node)
   ```
   Node i: [UX, UY, UZ, RX, RY, RZ] → global indices [6i, 6i+1, 6i+2, 6i+3, 6i+4, 6i+5]
   ```

3. **Boundary Condition Enforcement**
   - Identify fixed DOFs from supports
   - Extract partition: `K_ff` (free DOF stiffness), `F_f` (free DOF loads)
   - Solve reduced system, set fixed DOFs to zero displacement

---

## Element Types

### 1. **Beam Element** (3D Timoshenko/Euler-Bernoulli)
- **Degrees of Freedom**: 12 per element (6 per node × 2 nodes)
  - Translations: UX, UY, UZ
  - Rotations: RX, RY, RZ
- **Capabilities**: Full 3D bending, shear, torsion, axial
- **Applied with**: Member releases (partial fixity at ends)

### 2. **Truss/Cable Element** (Axial-only)
- **Degrees of Freedom**: 12 (6 per node, but only axial DOFs active)
- **Stiffness**: Only in local X direction
- **Special Feature**: Axial scale factor for cable slack analysis
  - `k[0,0] = (E×A/L) × axial_scale`
  - Reduces stiffness to model cables in compression

---

## Coordinate Transformation & Element Rotation

### Local-to-Global Transformation
**Transformation Matrix** (12×12 block diagonal):
```
T = [T3   0    0    0  ]
    [0    T3   0    0  ]
    [0    0    T3   0  ]
    [0    0    0    T3 ]
```
where `T3` is the 3×3 direction cosine matrix.

### Local Coordinate System Construction
For a member with start/end nodes and user-specified orientation vector:

1. **Local X**: Member axis
   ```
   x_local = (End - Start) / |End - Start|
   ```

2. **Local Y**: From user orientation
   ```
   y_proj = orientation - (orientation·x_local)×x_local
   y_local = y_proj / |y_proj|
   ```
   - User provides orientation vector (e.g., [0,1,0] for typical orientation)
   - Automatically orthogonalized to member axis

3. **Local Z**: Via cross product
   ```
   z_local = x_local × y_local
   ```

### Global Stiffness Matrix
```
K_global = T^T × k_local × T
```

---

## Stiffness Matrix Formulation

### Local Stiffness Matrix (12×12)
Standard **3D beam element** with uncoupled degrees of freedom:

#### Axial (DOF 0, 6)
$$k_{axial} = \frac{EA}{L}$$
```
[  EA/L              -EA/L           ]
[ -EA/L               EA/L           ]
```

#### Torsion (DOF 3, 9)
$$k_{torsion} = \frac{GJ}{L}$$
```
[  GJ/L              -GJ/L           ]
[ -GJ/L               GJ/L           ]
```

#### Bending about Z-axis (DOF 1, 5, 7, 11) — Y-displacements
**Hermite cubic shape functions** with 4 DOFs per axis:
$$k_{bending,z} = \frac{EI_z}{L^3} \begin{bmatrix}
12 & 6L & -12 & 6L \\
6L & 4L^2 & -6L & 2L^2 \\
-12 & -6L & 12 & -6L \\
6L & 2L^2 & -6L & 4L^2
\end{bmatrix}$$

#### Bending about Y-axis (DOF 2, 4, 8, 10) — Z-displacements
$$k_{bending,y} = \frac{EI_y}{L^3} \begin{bmatrix}
12 & -6L & -12 & -6L \\
-6L & 4L^2 & 6L & 2L^2 \\
-12 & 6L & 12 & 6L \\
-6L & 2L^2 & 6L & 4L^2
\end{bmatrix}$$

### Properties Used
- `E`: Young's modulus
- `G`: Shear modulus
- `A`: Cross-sectional area
- `Iy, Iz`: Second moments of inertia about y, z axes
- `J`: Polar moment of inertia (torsional constant)
- `L`: Member length

### Assumptions
- **Linear elastic material** (Hooke's law)
- **Small displacements** (linear kinematics)
- **No warping** (torsion via J only)
- **Decoupled DOFs** (axial, torsion, and two bending directions independent)

---

## Member End Releases

**Purpose**: Model partial fixity (pinned connections, etc.)

**Implementation**:
- 12-character string: "UXUYUZRXRYRZ" at start node + "UXUYUZRXRYRZ" at end node
- '0' = rigid connection, '1' = released (free)
- When released: zero out corresponding row/column, set diagonal to 1e-6 to prevent singularity

**Example**: `"000100000000"` releases Uy at start node.

---

## 1D Beam Analysis (Simplified)

For **1D beam problems** ([Load case solving in `beam1d/analysis.py`](src/beamy/beam1d/analysis.py)), a reduced FEM is used:

### Axial & Torsion
- **Element stiffness**: 2×2 linear springs
  $$k = \frac{EA}{L}\begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$$

### Transverse (Y, Z) Bending
- **Element stiffness**: 4×4 Hermite cubic (see above)
- **Shape functions**: Standard Hermite basis functions
  - Ensures $C^1$ continuity (displacement and slope continuous)
  - Exact for distributed loads on each element

---

## Analysis Types

### Frame Analysis
- Global FEM solve for all DOFs simultaneously
- Supports: 3D members, releases, different element types
- Output: Global displacements, member forces/moments via local stiffness

### 1D Beam Analysis
- Separate solves per direction (axial, torsional, bending-Y, bending-Z)
- Supported elements placed at problem-specific nodes
- Direct interpolation via Hermite basis for deflection profiles

---

## Matrix Properties

| Property | Value |
|----------|-------|
| Type | Symmetric |
| Positive Definite | Yes (before BC) |
| Sparse | Generally dense for frames |
| Bandwidth | Depends on node numbering |
| Condition Number | Improves with unit-consistent inputs |

---

## Common Method Classifications

**This is standard FEM using:**
- ✅ **Displacement Method** (stiffness method)
- ✅ **Direct Assembly** (no sub-structuring)
- ✅ **Euler-Bernoulli/Timoshenko Beam Theory**
- ✅ **Hermite Cubic Shape Functions** (transverse)
- ✅ **Linear Elements** (elastic)
- ✅ **Small Displacement Theory**
- ❌ No Galerkin/Variational (implicit in stiffness derivation)
- ❌ No iterative solvers (direct LU)
- ❌ No nonlinear effects

