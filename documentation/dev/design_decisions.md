# Design Decisions

This document tracks key architectural and design decisions made during the development of beamy.

## Coordinate System Reference Axis

**Decision:** The beam reference axis `(0,0)` corresponds to the **Centroid** of the cross-section, not the Shear Center.

### Context
In 1D beam analysis, the "axis" represents the line along which elements are defined. For symmetric sections (I-beams, RHS), the Centroid and Shear Center coincide. For asymmetric sections (Channels, Angles), they are offset.

We had to choose whether the reference axis `(0,0)` represents the Centroid or the Shear Center.

### Justification
We chose the **Centroid** as the reference axis for the following reasons:

1.  **Stress Calculations**: Fundamental stress formulas (`σ = N/A` and `σ = My/I`) are defined relative to the centroidal axis (neutral axis). Using the Shear Center would require transforming coordinates for every stress calculation.
2.  **Section Properties**: Standard libraries (`sectiony`) and handbooks define moments of inertia ($I_y, I_z$) about the centroid. Using the Shear Center would require applying the Parallel Axis Theorem for every element property.
3.  **Axial Loading Intuition**: In structural engineering, axial members are typically modeled as being loaded through their geometric center (centroid) to avoid induced moments.

### Implications & Handling
*   **Transverse Loading**: Applying a shear force at the Centroid of an asymmetric section physically induces torsion.
    *   *Solution*: Beamy's solver (`analysis.py`) automatically calculates and applies this induced torsional moment ($T = -SC_y F_z + SC_z F_y$) using the section's Shear Center offset properties.
*   **Visual Inconsistency**: The 3D plotter shifts the section so the Shear Center aligns with the visual axis to simplify twist visualization. This is a visual-only transformation and does not affect the underlying physics calculation.

