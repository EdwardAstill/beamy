Below are **theoretical structural mechanics notes on Euler–Bernoulli beam theory**, written at the level typically used for continuum-to-beam reduction and FEM formulation.

---

## Euler–Bernoulli Beam Theory

Euler–Bernoulli beam theory is a classical one-dimensional structural model used to describe the bending behaviour of slender members. It is derived from linear elasticity under restrictive kinematic assumptions that eliminate shear deformation effects.

---

### 1. Scope and Applicability

Euler–Bernoulli theory is valid when:

* Beam length ( L ) is much larger than cross-sectional dimensions.
* Bending deformation dominates.
* Transverse shear deformation is negligible.
* Material behaviour is linear elastic.
* Rotations and strains are small.

Typical applications include slender steel or concrete beams under service-level loading.

---

### 2. Kinematic Assumptions

The defining assumption of Euler–Bernoulli theory is:

> **Plane cross-sections remain plane and normal to the deformed beam axis.**

Consequences:

* No relative sliding between material layers.
* Shear strain ( \gamma_{xz} = 0 ).
* Cross-section rotation equals the slope of the deflection curve.

For transverse displacement ( w(x) ):
[
\theta(x) = \frac{dw}{dx}
]

---

### 3. Displacement Field

For a beam aligned with the (x)-axis and bending in the (x\text{–}z) plane:

[
u_x(x,z) = -z \frac{dw}{dx}
]
[
u_z(x) = w(x)
]

There is no independent rotation field.

---

### 4. Strain–Displacement Relations

The only nonzero strain component is axial strain:

[
\epsilon_{xx} = \frac{\partial u_x}{\partial x}
= -z \frac{d^2 w}{dx^2}
]

This introduces curvature:
[
\kappa = \frac{d^2 w}{dx^2}
]

---

### 5. Constitutive Law

Linear elastic material behaviour:
[
\sigma_{xx} = E \epsilon_{xx}
]

Stress varies linearly through the depth of the beam.

---

### 6. Stress Resultants

#### 6.1 Bending Moment

[
M(x) = \int_A \sigma_{xx} z , dA
= -E \kappa \int_A z^2 , dA
= EI \kappa
]

where:
[
I = \int_A z^2 , dA
]

#### 6.2 Shear Force

Shear force is obtained from equilibrium, not constitutive relations:
[
V(x) = \frac{dM}{dx}
]

Shear stress distribution is not explicitly modeled.

---

### 7. Governing Differential Equation

From equilibrium of an infinitesimal beam element:

[
\frac{d^2 M}{dx^2} = q(x)
]

Substituting ( M = EI \kappa ):

[
EI \frac{d^4 w}{dx^4} = q(x)
]

This is the **Euler–Bernoulli beam equation**.

---

### 8. Boundary Conditions

Each beam end provides two boundary conditions, chosen from:

* Deflection ( w )
* Rotation ( \theta = w' )
* Bending moment ( M = EI w'' )
* Shear force ( V = EI w''' )

Examples:

* Clamped end: ( w = 0, ; w' = 0 )
* Simply supported end: ( w = 0, ; M = 0 )

---

### 9. Energy Interpretation

The strain energy of bending is:

[
U = \frac{1}{2} \int_0^L EI \left( \frac{d^2 w}{dx^2} \right)^2 dx
]

Euler–Bernoulli theory is variationally consistent and forms the basis of beam finite elements via the principle of minimum potential energy.

---

### 10. Relation to Finite Element Formulation

* Leads to ( C^1 )-continuous displacement fields.
* Classical beam elements use cubic Hermite shape functions.
* Degrees of freedom per node:

  * Transverse displacement ( w )
  * Rotation ( \theta )

---

### 11. Limitations

Euler–Bernoulli theory:

* Neglects transverse shear deformation.
* Overestimates stiffness for deep or short beams.
* Is inaccurate for sandwich or low-shear-modulus materials.
* Cannot capture shear locking effects because shear is absent.

For such cases, **Timoshenko beam theory** is required.

---

### Summary

Euler–Bernoulli beam theory models bending by assuming shear-free kinematics and linear elastic material behaviour. It reduces a 3D continuum problem to a 1D formulation with high accuracy for slender beams, forming the theoretical foundation of classical beam analysis and many finite element implementations.
