
# Analysis (post-upgrade: second-order + imperfections + analysis modes)

This document describes how **beamy**’s *frame* analysis is intended to work once the planned stability upgrades are implemented:

- a common **second-order elastic** analysis engine (P–Δ + P–δ),
- a single, consistent **equivalent imperfection** mechanism (typically notional loads),
- **code-specific analysis modes** (AISC 360 DAM, Eurocode 3 global imperfections, AS 4100 second-order / amplification triggers),
- optional **6-DOF nodal springs** to model semi-rigid restraint.

The overall philosophy is:

> Put stability physics in the analysis, keep design checks as capacity comparisons.

That lets member checks (AISC Chapter E/F/H, future EC3/AS checks) consume demands that already include the relevant second-order amplification.

---

## Table of contents

- [What “analysis” means in beamy](#what-analysis-means-in-beamy)
- [Analysis modes](#analysis-modes)
- [Core engine (elastic second-order)](#core-engine-elastic-second-order)
- [Imperfections (equivalent)](#imperfections-equivalent)
- [Stiffness rules (code wrappers)](#stiffness-rules-code-wrappers)
- [Boundary conditions, releases, and numerical stabilization](#boundary-conditions-releases-and-numerical-stabilization)
- [Cables and other nonlinear element behaviors](#cables-and-other-nonlinear-element-behaviors)
- [Node springs (6-DOF)](#node-springs-6-dof)
- [Results and downstream checks](#results-and-downstream-checks)
- [Limitations](#limitations)

---

## What “analysis” means in beamy

beamy has two closely related analysis layers:

1. **Single-member analysis** via `LoadedMember` (implemented as a 2-node / 1-member frame solve).
2. **Frame** analysis (3D frame), used to solve connected structures with 6 DOFs per node.

The current implementation prioritizes **direct frame analysis** as the canonical approach:

**Key principle:**

> Frame analysis is the authoritative source for member demands. Member design checks consume those demands directly via equilibrium without re-solving members as isolated beams.

This ensures:
- **Consistency**: all checks use the same demand source (frame-recovered end forces + distributed loads)
- **Correctness**: symmetric members yield symmetric demands (no artificial boundary condition asymmetries)
- **Efficiency**: no redundant 1D FEM solves; demands computed via equilibrium

### Demand computation path (Direct Frame Analysis)

For each member, beamy:
1. Recovers end forces from the frame analysis (local axes)
2. Collects distributed loads applied to that member (converted to local coordinates if needed)
3. Computes internal actions (N, V, M) at arbitrary points using equilibrium without re-solving
4. Assembles a `MemberActionProfile` for use in design checks

This approach automatically handles:
- Members split at intermediate nodes (stitches segments seamlessly)
- Linearly varying distributed loads
- Asymmetric member orientations and end conditions
- Symmetric structures with symmetric results

---

## Analysis modes

beamy will expose “analysis modes” as explicit, user-visible choices. Each mode maps to a consistent recipe:

- **GENERIC_LINEAR**
	- First-order (small displacement) linear elastic.
	- No P–Δ/P–δ.
	- No required imperfections.

- **SECOND_ORDER_ELASTIC**
	- Second-order elastic analysis via geometric stiffness and iteration.
	- User chooses imperfection model (typically notional loads).
	- User chooses stiffness modifiers (if any).

- **AISC360_DAM** (Direct Analysis Method)
	- Second-order elastic analysis + imperfections (usually notional loads).
	- AISC-specific stiffness reduction rules (e.g., reduced bending stiffness).
	- When converged and conditions are met, compression checks can default to **K = 1.0**.

- **EC3_GLOBAL_IMPERFECTIONS**
	- Global analysis includes imperfections when relevant.
	- Second-order may be required depending on stability criteria or user override.
	- Imperfection magnitudes/format follow EC3 conventions.

- **AS4100_SECOND_ORDER**
	- Second-order elastic analysis for cases that require it.
	- Otherwise, an amplified-first-order path may be used (future enhancement).
	- Notional horizontal force convention is represented via the notional-load model.

Each analysis run stores its selected mode/settings on the results so that checkers can enforce consistency.

---

## Core engine (elastic second-order)

### Summary

The engine solves equilibrium with geometric nonlinearity effects using an iterative sequence:

1. Assemble material stiffness $K$ (linear elastic, small strain).
2. Solve $K u = F$.
3. Recover element axial forces $P$ (and member end forces) from the displacement solution.
4. Assemble geometric stiffness $K_g(P)$.
5. Solve $(K + K_g) u = F$.
6. Iterate until convergence.

This captures:

- **P–Δ** (global sway amplification)
- **P–δ** (member curvature amplification)

### Sign conventions (critical for correctness)

Second-order behavior depends on the sign of axial force in the geometric stiffness formulation.

Requirement:

- **Compression** must be destabilizing (reduces effective lateral stiffness).
- **Tension** must be stabilizing.

The implementation will define a single sign convention for member axial force recovery (local axes) and use that same convention when building $K_g(P)$.


### Convergence behavior

The solve is iterative and must converge to be considered “checker-safe”. The results will record:

- `converged: bool`
- `iterations: int`
- one or more convergence metrics (e.g., relative displacement change)

If the second-order solve does **not** converge, beamy will raise a clear error and report the last iterate for debugging.

### Convergence hardening (load stepping + damping)

Near-instability frames are exactly where second-order iteration is most valuable and also most likely to struggle.

To improve robustness, the solver will optionally support:

- **Load stepping**: apply the load case in $n$ increments and converge each increment.
- **Under-relaxation**: damp oscillatory updates using a relaxation factor $\omega$.

These features are optional, but when used they will be recorded in the analysis results metadata.

---

## Imperfections (equivalent)

The analysis supports imperfections using a single consistent mechanism.

### Notional loads (recommended default)

Notional loads are implemented as equivalent lateral nodal forces derived from vertical design load:

$$H = \alpha \; P$$

where:

- $\alpha$ is the notional-load factor (often $0.002$),
- $P$ is a gravity measure (per node or per “level”, depending on the model).

Implementation intent:

- A preliminary gravity-only solve (linear is acceptable) provides reactions / vertical force measures.
- Lateral notional forces are generated from that measure and added to the active load case.

### Deterministic definition of P (for notional loads)

To avoid sensitivity to support stiffness (especially once springs exist), the default notional-load generator will use an **input-based gravity measure**:

- applied vertical nodal loads, plus
- member distributed vertical loads lumped to nodes or levels using a documented scheme.

If a reaction-based method is used (optional), it must be explicitly selected and recorded in results metadata.

### Initial sway geometry (optional / future)

As an alternative to notional loads, an explicit out-of-plumb geometry can be supported by applying small node offsets before analysis.

---

## Stiffness rules (code wrappers)

The common engine accepts **stiffness modifiers** that can be applied per element.

Examples of modifier knobs:

- bending stiffness multiplier (e.g., reduced $EI$)
- optional axial stiffness multiplier ($EA$) if required by a specific methodology

Wrappers (AISC/EC3/AS) choose the modifiers and imperfection settings, but the solver path stays the same.

---

## Boundary conditions, releases, and numerical stabilization

### Supports and constraints

- Nodes have 6 DOFs (3 translations + 3 rotations).
- Nodal supports constrain DOFs using boolean fixity.
- Member end constraints (if present) also constrain DOFs at member ends.
- Member releases represent end releases in the element formulation (local stiffness terms removed / modified).

### Automatic stabilization (important)

Frame models may legitimately contain nodes connected only by axial elements (truss/cable). Rotations at those nodes can have no stiffness contribution, which can make the global stiffness matrix singular.

To keep the system solvable, beamy uses a **minimal numerical stabilizer**:

- rotations (RX/RY/RZ) at truss/cable-only nodes are automatically fixed.

This is a numerical device, not a physical restraint.

---

## Cables and other nonlinear element behaviors

Tension-only cables are handled using an outer iteration:

1. Solve for displacements.
2. Recover cable axial force.
3. If a cable is in compression, reduce its axial stiffness to a near-zero “slack” value.
4. Repeat until the active set stops changing.

When second-order analysis is enabled, the intended structure is:

- **Outer loop**: cable slack/taut set changes
- **Inner loop**: second-order convergence for a fixed cable state

This keeps the nonlinearities controlled and makes convergence behavior easier to diagnose.

---

## Node springs (6-DOF)

To model semi-rigid restraint (partial fixity) without forcing pinned/fixed assumptions, beamy will support uncoupled nodal springs:

- translational: $k_x, k_y, k_z$
- rotational: $k_{rx}, k_{ry}, k_{rz}$

Implementation intent:

- springs contribute directly to the global stiffness matrix (typically diagonal terms in the uncoupled case),
- compatible with both first-order and second-order analysis.

This feature improves realism for *all* standards because joint restraint strongly influences sway and end moments.

### Future-proofing: coupled springs

Uncoupled (diagonal) springs are the initial implementation.
The API will be designed so a node can eventually accept a full **6×6 coupled spring matrix** (even if only diagonal terms are supported at first).

---

## Results and downstream checks

### Stored results

Frame analysis will provide (at minimum):

- nodal displacements per node
- reactions at constrained nodes
- member end forces per member (local)
- convenience APIs to obtain member internal actions along the member (using equilibrium from recovered end forces + distributed loads)

### Design-grade status

Beyond `converged`, the analysis results will include a stricter design-consumption flag:

- `design_grade_ok: bool`
- `design_grade_notes: list[str]`

Examples of reasons `design_grade_ok` might be false:

- nonconvergence
- load stepping reached maximum without convergence
- numerical stabilization DOFs were added (still solvable, but not “pure physics”)


### Results metadata (critical for checkers)

Each analysis run should store:

- `analysis_method` (e.g., `AISC360_DAM`)
- whether second-order was enabled
- imperfection model and factor
- stiffness rule set applied
- convergence information

Design checkers can then apply method-specific assumptions safely. Example:

- AISC Chapter E can default to $K = 1.0$ only when analysis method is DAM and the solve converged.

### Canonical member demand API (for thin checkers)

To keep AISC/EC3/AS wrappers thin and consistent, beamy will provide one canonical demand interface per member, using **local axes**:

- end forces at both ends: $(N, V_y, V_z, T, M_y, M_z)$
- internal action stations or envelopes along the member for the same quantities
- one documented sign convention referenced across all checkers

---

## Limitations

This analysis approach is intentionally “elastic + stability-focused”, not a full nonlinear FEM program.

Current/expected scope limits:

- elastic material behavior (no plastic hinges)
- small strain; second-order captured via geometric stiffness iteration
- no full large-rotation corotational beam formulation (initially)
- imperfections represented by notional loads and/or prescribed initial sway (not fabrication-residual stress modeling)

These choices are deliberate: they target repeatable, checker-friendly demands for steel design workflows.

