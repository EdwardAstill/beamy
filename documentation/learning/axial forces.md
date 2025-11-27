The math behind `solve_axial_forces`

## Statically determinate cases (one axial support)

If there is only one axial restraint, the axial force distribution is statically determinate.

- First, find the reaction force at the axial support from global equilibrium.

- Then move along the beam from the free end towards the support. At each station, if you cross a point force (including forces at supports), update the internal normal force:

  - A positive external axial force reduces the internal normal force.

  - A negative external axial force increases the internal normal force.

At a truly free axial end, the internal normal force is zero.

## Statically indeterminate cases (more than one axial support)

If there is more than one axial restraint, the problem is statically indeterminate in the axial direction and we must use the stiffness method.

The idea is to set up simultaneous equations so that:

- Equilibrium is satisfied at each node.

- Compatibility is satisfied through the element stiffness relationships.

- Support conditions (fixed/free/settlement) are enforced via boundary conditions.

There are two types of nodes in the stiffness method:

- **Free nodes**

  - Unknown displacements

  - Known nodal forces (including zero if no external load)

- **Support nodes**

  - Known displacements (usually 0)

  - Unknown forces (reactions)

### Element stiffness matrix

Between every pair of nodes there is a 1D bar (element) with stiffness

$$k_e = \frac{E A}{L_e}$$

For an element between node $i$ and node $j$:

- Nodal displacements: $u_i, u_j$

- Nodal forces: $F_i, F_j$

The element equations are:

$$\begin{aligned}
F_i &= k_e (u_i - u_j) \\
F_j &= k_e (u_j - u_i)
\end{aligned}$$

In matrix form:

$$\begin{bmatrix}
F_i \\
F_j
\end{bmatrix}
=
k_e
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}
\begin{bmatrix}
u_i \\
u_j
\end{bmatrix}$$

This can also be derived using the standard FE procedure:

$$u(x) = N_1 u_i + N_2 u_j, \quad
\varepsilon = B u, \quad
K_e = \int B^T E A B \, dx$$

### Global stiffness matrix

We assemble the global stiffness matrix by superimposing the contributions from each element into the correct rows and columns corresponding to the global node numbers.

Consider 3 nodes (1, 2, 3) connected by 2 elements with stiffnesses $k_1$ (between 1–2) and $k_2$ (between 2–3). Let the nodal displacements be $u_1, u_2, u_3$ and nodal forces $F_1, F_2, F_3$. The equilibrium equations at each node are:

$$\begin{aligned}
F_1 &= k_1 (u_1 - u_2) \\
F_2 &= k_1 (u_2 - u_1) + k_2 (u_2 - u_3) \\
F_3 &= k_2 (u_3 - u_2)
\end{aligned}$$

In matrix form:

$$\begin{bmatrix}
F_1 \\
F_2 \\
F_3
\end{bmatrix}
=
\begin{bmatrix}
k_1 & -k_1 & 0 \\
-k_1 & k_1 + k_2 & -k_2 \\
0 & -k_2 & k_2
\end{bmatrix}
\begin{bmatrix}
u_1 \\
u_2 \\
u_3
\end{bmatrix}$$

For $n$ nodes with element stiffnesses $k_1, \dots, k_{n-1}$, the global stiffness matrix $K$ is tridiagonal. Each node only interacts with its immediate neighbours, so the non-zero entries occur on the main diagonal and the two adjacent diagonals.

### Apply boundary conditions

We now apply the support and loading conditions:

- At **support nodes**: the displacement $u$ is known (typically 0), and the reaction force is unknown.

- At **free nodes**: the nodal axial force is known (external load, or zero if unloaded), and the displacement is unknown.

We reorder the system so that the unknown displacements $u_f$ (free DOFs) come first, and the known displacements $u_s$ (support DOFs) come last:

$$\begin{bmatrix}
F_f \\
F_s
\end{bmatrix}
=
\begin{bmatrix}
K_{ff} & K_{fs} \\
K_{sf} & K_{ss}
\end{bmatrix}
\begin{bmatrix}
u_f \\
u_s
\end{bmatrix}$$

Here:

- $u_s$ is known from the support conditions.

- $F_f$ is known from the applied loads (including zeros at unloaded free nodes).

We solve

$$K_{ff} u_f = F_f - K_{fs} u_s$$

for the unknown displacements $u_f$. Then we back-substitute to find the reaction forces at the supports:

$$F_s = K_{sf} u_f + K_{ss} u_s$$

Once the nodal forces are known, the internal axial force in each element can be recovered from the element equations.
