
## ---

**1\. Structural Mechanics (The Engineering Physics)**

This category defines the "Rules of the Universe" for your model. It tells the software how a member should behave under load before any math is applied.

* **Elasticity & Constitutive Laws**: Understanding the relationship between stress and strain (Hooke’s Law) for 1D elements. In Beamy, this translates to how you use **Young’s Modulus ($E$)** and **Shear Modulus ($G$)** to define stiffness.  
* **Beam Theories**:  
  * **Euler-Bernoulli**: For slender members where bending dominates.  
  * **Timoshenko-Ehrenfest**: Necessary if you want to account for shear deformation in deep beams.  
* **Failure & Stability Theory**:  
  * **Yielding**: When the material itself fails.  
  * **Buckling**: When the *geometry* fails due to instability (crucial for your design/ module).  
* **Geometric Nonlinearity**: The theory behind **P-Delta** effects, where the displaced shape of the structure creates additional internal moments.

## ---

**2\. Finite Element Method (The Numerical Engine)**

FEM is the bridge that turns the physics of "continuous" beams into "discrete" computer-readable matrices.

* **The Direct Stiffness Method**: The primary algorithm for frame analysis. You must understand how to build the $12 \\times 12$ **Element Stiffness Matrix** for a 3D beam.  
* **Coordinate Transformations**: Using rotation matrices to convert local member coordinates (along the beam) to global frame coordinates (X, Y, Z).  
* **Static Condensation**: The specific numerical trick needed for your releases feature. It allows you to "zero out" stiffness at pinned ends while maintaining the integrity of the global matrix.  
* **Numerical Linear Algebra (The Solver)**:  
  * **Linear**: Solving $Ku \= F$ using LU Decomposition or Sparse solvers.  
  * **Nonlinear**: Implementing the **Newton-Raphson** iteration to solve second-order cases by updating the "Tangent Stiffness" matrix until convergence.

## ---

**3\. Engineering Design (The Compliance Layer)**

This is where you apply human-defined safety standards to the raw physics results.

* **Load Factoring (LRFD/ASD)**: Understanding how to scale loads and strengths based on statistical safety (e.g., $1.2D \+ 1.6L$).  
* **Member Capacity Checks**:  
  * **AISC 360 (Steel)**: Calculating $\\phi P\_n$ (axial strength) and $\\phi M\_n$ (bending strength).  
  * **Interaction Equations**: The logic that checks if a beam can handle *both* bending and axial force simultaneously without failing.  
* **Limit States**:  
  * **Strength**: Does it break?.  
  * **Serviceability**: Does it vibrate or deflect too much for the user's comfort?.  
* **Unbraced Length Theory**: The mechanics of **Lateral-Torsional Buckling**, which dictates how much bending a beam can take before it twists out of plane.

### ---

**How they interact in your code**

| Feature | Structural Mechanics | FEM | Engineering Design |
| :---- | :---- | :---- | :---- |
| **Pinned End** | Equilibrium at joint. | Static Condensation. | End-connection detailing. |
| **Distributed Load** | Fixed-end moments. | Equivalent Nodal Loads. | Load path validation. |
| **Utilization Ratio** | Stress calculation. | Displacement recovery. | Capacity vs. Demand ratio. |

**Would you like me to focus on the Python code for the $12 \\times 12$ matrix assembly (FEM), or the specific AISC formulas for steel capacity (Design)?**