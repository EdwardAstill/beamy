# Design Module Reference

The `design/` module performs code checks on analysis results. It is independent of solver internals.

## `design/aisc/`

Public API for AISC design checks (e.g., AISC 360).

### `design/aisc/compression.py`

**Chapter E** (Compression Members).

*   Computes nominal compressive strength $P_n$.
*   Calculates design strength $\phi P_n$ (LRFD) or allowable strength $P_n / \Omega$ (ASD).
*   Requires:
    *   Section and material properties.
    *   Effective lengths ($K L$).

### `design/aisc/flexure.py`

**Chapter F** (Flexural Members).

*   Computes nominal flexural strength $M_n$ for strong and weak axes.
*   Handles **Lateral-Torsional Buckling (LTB)** checks.
*   Requires:
    *   Unbraced length $L_b$.
    *   Moment gradient factor $C_b$.

### `design/aisc/interaction.py`

**Chapter H** (Interaction).

*   Combines axial ($P_u$) and flexural ($M_{ux}, M_{uy}$) demands.
*   Uses capacities calculated in Chapters E and F.
*   Returns:
    *   Utilization ratio.
    *   Governing equation reference.

## Design Workflow

Design functions accept:
1.  **Demands**: From `MemberResult` or `FrameResult`.
2.  **Design Parameters**: User-provided values (K factors, unbraced lengths, Cb) stored in `Member.design` properties.
    *   **Validation**: These properties should be pre-validated by `Frame.validate()` to ensure defaults (e.g., $L_b = L$) are set if omitted.

