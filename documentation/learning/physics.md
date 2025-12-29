# Physics and Analysis Methods

This document outlines the theoretical foundations and numerical methods used in beamy for structural beam analysis.

## Overview

Beamy uses the **Finite Element Method (FEM)** to solve for displacements and reactions in 1D beam structures. The analysis is split into four independent problems:

1. **Axial** (tension/compression along x-axis)
2. **Torsion** (twisting about x-axis)
3. **Transverse Bending in Y** (loads in y-direction, bending about z-axis)
4. **Transverse Bending in Z** (loads in z-direction, bending about y-axis)

Each problem is solved separately using appropriate element formulations.

---

## Beam Theory: Euler-Bernoulli

Transverse bending analysis uses **Euler-Bernoulli beam theory**, which assumes:

- Plane sections remain plane after deformation
- Plane sections remain perpendicular to the neutral axis (no shear deformation)
- Small deflections and rotations
- Linear elastic material behavior

This theory is valid for slender beams where the length is much greater than the cross-sectional dimensions (L/d > 10 typically).

The governing equation is:

```
EI * d⁴w/dx⁴ = q(x)
```

where:
- `E` = Young's modulus
- `I` = Second moment of area about the bending axis
- `w` = Transverse deflection
- `q(x)` = Distributed load

---

## Finite Element Method

### Element Types

**Axial and Torsion Elements (1 DOF per node)**

Uses linear (2-node) bar elements with shape functions:
```
N₁ = 1 - ξ
N₂ = ξ
```
where ξ ∈ [0, 1] is the normalized coordinate.

Local stiffness matrix:
```
k = (1/L) * [[ 1, -1],
             [-1,  1]]
```

**Bending Elements (2 DOFs per node: deflection w and rotation θ)**

Uses **Hermite cubic** beam elements, which provide C¹ continuity (continuous displacement and slope).

Hermite shape functions:
```
N₁ = 1 - 3ξ² + 2ξ³       (deflection at node 1)
N₂ = L(ξ - 2ξ² + ξ³)     (rotation at node 1)
N₃ = 3ξ² - 2ξ³           (deflection at node 2)
N₄ = L(-ξ² + ξ³)         (rotation at node 2)
```

Local stiffness matrix (4×4):
```
k = (EI/L³) * [[ 12,    6L,   -12,    6L  ],
               [  6L,   4L²,  -6L,   2L² ],
               [-12,   -6L,    12,   -6L  ],
               [  6L,   2L²,  -6L,   4L² ]]
```

### Global Assembly

Element stiffness matrices are assembled into a global stiffness matrix K using standard FEM assembly:

```
K[global_dofs] += k[local_dofs]
```

### Boundary Conditions

Supports are defined using a 6-digit string representing constraints on 6 DOFs:
- Position 1: Ux (axial translation)
- Position 2: Uy (transverse y)
- Position 3: Uz (transverse z)
- Position 4: Rx (torsional rotation)
- Position 5: Ry (rotation about y)
- Position 6: Rz (rotation about z)

Where `1` = fixed, `0` = free.

Boundary conditions are applied by partitioning the system:
```
[K_ff  K_fc] [d_f]   [f_f]
[K_cf  K_cc] [d_c] = [f_c]
```

where subscript `f` = free DOFs, `c` = constrained DOFs.

Solving: `d_f = K_ff⁻¹ * f_f`

### Reaction Calculation

Reactions are computed as:
```
r = K * d - f
```

where `r` at constrained DOFs gives the support reactions.

---

## Automatic Node Insertion

To ensure accurate deflection calculations, beamy automatically inserts internal "free" nodes at load application points. This is necessary because the FEM computes displacements only at nodes.

For a simply supported beam with a mid-span load, without an internal node the solver would only compute displacements at the supports (which are zero by definition). By inserting a free node at the load location, the solver can compute the actual deflection there.

---

## Internal Force Diagrams

Shear force and bending moment diagrams are computed by **load accumulation**:

**Shear Force V(x)**:
```
V(x) = Σ F_i  for all loads at positions x_i ≤ x
```

**Bending Moment M(x)**:
```
M(x) = Σ F_i * (x - x_i) + Σ M_j  for all loads/moments at positions ≤ x
```

This is equivalent to integrating the load diagram from the left support.

---

## Displacement Interpolation

### Bending Deflection

Deflections between nodes are interpolated using **Hermite shape functions**:

```
w(ξ) = N₁*w₁ + N₂*θ₁ + N₃*w₂ + N₄*θ₂
```

This ensures smooth, physically realistic deflection curves that match the cubic polynomial assumed in Euler-Bernoulli theory.

### Axial/Torsional Displacement

Uses **linear interpolation** between nodal values:

```
u(x) = numpy.interp(x, x_nodes, d_nodes)
```

---

## Stress Calculations

### Axial Stress
```
σ = N / A
```
where N = axial force, A = cross-sectional area.

### Bending Stress
```
σ = M * c / I
```
where M = bending moment, c = distance to extreme fiber, I = second moment of area.

### Shear Stress (Average)
```
τ = V / A
```
where V = shear force. This is a simplified average; actual shear stress varies parabolically across the section.

### Torsional Shear Stress
```
τ = T * r / J
```
where T = torque, r = radial distance to extreme fiber, J = polar moment of area.

### Von Mises Stress

A conservative estimate combining all stress components:

```
σ_vm = √(σ_total² + 3 * τ_total²)
```

where:
- `σ_total = |σ_axial| + |σ_bending_y| + |σ_bending_z|`
- `τ_total = |τ_shear_y| + |τ_shear_z| + |τ_torsion|`

This is conservative because it assumes all maximum stresses occur at the same point, which may not be physically accurate for all cross-sections.

---

## Coordinate System

Beamy uses a right-handed coordinate system where the **reference axis (0,0) corresponds to the Centroid of the cross-section**.

- **x**: Along the beam axis (longitudinal)
- **y**: Vertical (typically)
- **z**: Horizontal (typically)

Bending in the x-y plane produces moments about z (Mz).
Bending in the x-z plane produces moments about y (My).

*(See `documentation/dev/design_decisions.md` for more details on this choice)*

---

## Limitations

1. **Euler-Bernoulli assumption**: Shear deformation is neglected. For deep beams (L/d < 10), Timoshenko beam theory would be more accurate.

2. **Linear elastic**: Material is assumed to behave linearly. No plasticity or yielding.

3. **Small deflections**: Large deflection (geometric nonlinearity) effects are not considered.

4. **Prismatic beams**: Cross-section is constant along the length.

5. **Static analysis**: Dynamic effects (vibration, impact) are not included.

---

## References

- Euler-Bernoulli Beam Theory
- Finite Element Method for Structural Analysis
- Hermite Interpolation Polynomials
- Matrix Structural Analysis

