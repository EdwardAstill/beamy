
# Common second-order analysis engine (AISC + EC3 + AS 4100) — implementation outline

Goal: implement one reusable **second-order elastic frame analysis** engine (P–Δ + P–δ) with a consistent **equivalent imperfection mechanism** (notional loads or initial sway), then layer **code-specific stiffness + rules** on top via “analysis modes”.

This is intended to be the common baseline for:

- **AISC 360**: Direct Analysis Method (DAM)
- **Eurocode 3 (EN 1993-1-1)**: global analysis including imperfections; 1st vs 2nd order triggers
- **AS 4100**: second-order vs amplified-first-order triggers; notional horizontal force convention

---

## 1) Architecture (separate engine from wrappers)

### 1.1 Add an explicit analysis-mode concept

Add an enum/string to represent the chosen analysis track:

- `GENERIC_LINEAR` (current)
- `ELASTIC_SECOND_ORDER` (common engine with user-selected imperfection + stiffness rules)
- `AISC360_DAM`
- `EC3_GLOBAL_IMPERFECTIONS`
- `AS4100_SECOND_ORDER`

Store it in the analysis results so checkers know what assumptions were used.

### 1.2 Define a single settings object

Create `FrameAnalysisSettings` (dataclass) to drive all analysis runs:

- `method`: one of the modes above
- `second_order`: bool (derived from method)
- `imperfection`: `"none" | "notional_loads" | "initial_sway"`
- `notional_factor`: default `0.002` (wrapper can override)
- `stiffness_rules`: `"none" | "aisc_dam" | "ec3" | "as4100"`
- `max_iter`, `tol_u_rel`, `tol_u_abs`
- `use_geometric_stiffness`: bool
- `include_p_delta`, `include_p_small_delta` (can be one combined knob initially)
- `cable_tension_only`: keep existing behavior

Deliverable: settings round-trip into `LoadedFrame` and appear in reports.

---

## 2) Core engine: second-order elastic analysis

### 2.1 Minimum viable second-order strategy

Implement **iterative second-order** using geometric stiffness:

1. Assemble material stiffness `K(u)` (initially constant, linear elastic)
2. Solve `K u = F` (existing)
3. Recover member axial forces `P` (existing force recovery)
4. Assemble geometric stiffness `Kg(P)` per element
5. Solve `(K + Kg) u = F`
6. Iterate until converged

Notes:

- Start with frame “beam” elements only; treat truss/cable as axial-only geometric stiffness (or omit initially).
- Keep the cable tension-only loop, but structure loops so behavior is deterministic:
	- Outer loop: cable active set changes
	- Inner loop: second-order convergence for fixed cable stiffness state

### 2.2 Where it likely lands in the current code

- The entry point is `LoadedFrame._analyze_frame()` in src/beamy/frame/analysis.py.
- Split out analysis into composable functions:
	- `_assemble_K(frame, node_to_idx, member_axial_scales, stiffness_modifiers) -> (K, per_member_mats)`
	- `_recover_member_end_forces(per_member_mats, displacements) -> member_end_forces`
	- `_assemble_Kg(frame, node_to_idx, member_end_forces, ...) -> Kg`

### 2.3 Geometric stiffness (Kg) implementation detail

Implement `Kg` at the element level:

- For beam-column: geometric stiffness based on current axial force `P` in the element.
- Transform `Kg_local -> Kg_global` using existing transformation matrix machinery.

Deliverable: `analyze_frame_geometry(...)` grows an optional `include_geometric_stiffness` path or a sibling function.

#### 2.3.1 Critical: tension vs compression sign conventions

Define and document a single sign convention for recovered axial force $P$ and ensure the geometric stiffness uses it correctly:

- **Compression** must be destabilizing (effective lateral stiffness decreases)
- **Tension** must be stabilizing

Implementation requirement:

- write one unit test or small benchmark that applies the same lateral load to (a) a column in compression and (b) a tie in tension, and confirms the expected trend in lateral stiffness.
- store a short note in results metadata that records the sign convention used (e.g., “local +x axial force sign”).

### 2.4 Convergence + failure behavior

Add diagnostics to `LoadedFrame`:

- `converged: bool`
- `iterations: int`
- `residual_norms: list[float]` (optional)
- `warnings: list[str]`

Stop criteria (pick one consistent rule early):

- relative displacement change: `||u_i - u_{i-1}|| / max(||u_i||, eps) < tol_u_rel`
and/or
- absolute displacement change: `||u_i - u_{i-1}|| < tol_u_abs`

If not converged: raise with a clear message and keep the last iterate for debugging.

#### 2.4.1 Add robustness: load stepping + damping

Add optional continuation controls (off by default) for hard sway cases:

- **Load stepping**: solve the same load case in $n$ increments (e.g., 5–20 steps), carrying forward $u$ as the initial guess.
- **Under-relaxation** (simple damping): if iteration oscillates, update $u_{k+1} \leftarrow (1-\omega)u_k + \omega u_{k+1}$ with $0 < \omega \le 1$.

Deliverables:

- `n_steps` setting
- `relaxation_omega` setting
- results store whether stepping/damping was used and whether the solve would have failed without it.

---

## 3) Equivalent imperfections (one consistent mechanism)

### Option A (recommended first): notional lateral loads

Implement a helper that creates notional loads from vertical design load:

- compute “gravity compression measure” per node or per level
- apply `H = notional_factor * P` in X and/or Y

Implementation approach (simple + practical):

1. Run a preliminary gravity-only solve (linear is OK)
2. Use reactions (or a vertical load sum) to get `P_level` or `P_node`
3. Apply lateral forces at those nodes

#### 3.A.1 Make notional loads deterministic: define P source

Choose and codify a single default source for $P$ so notional loads are reproducible and not sensitive to support stiffness modeling:

- Default recommendation: **input-based gravity measure** (from applied vertical nodal forces + member distributed vertical loads lumped to nodes/levels).

If a reaction-based measure is ever used (optional):

- it must be explicitly selected and recorded in results metadata.

Deliverables:

- `FrameLoadCase.with_notional_loads(settings) -> FrameLoadCase` (pure function)
- or `LoadedFrame._build_effective_load_case()` that returns `loads + notional(loads)`

### Option B (future): initial sway geometry

Implement a geometry perturbation mechanism:

- node position offsets before assembly
- a clear way to “reset to original geometry”

This is optional; if notional loads exist, most workflows can rely on them.

---

## 4) Code-specific stiffness rules (wrappers)

### 4.1 Stiffness modifiers should be a first-class input

Create a “stiffness modifier” concept that affects element matrices during assembly:

- `E_mod` multiplier
- `Iyy_mod`, `Izz_mod` multipliers (or direct `EI` multipliers)
- optional axial `EA_mod`

Apply it inside element stiffness creation in the solver path.

### 4.2 AISC360_DAM wrapper

Wrapper responsibilities:

- enable second-order
- enable imperfections (typically notional loads)
- apply DAM stiffness reduction policy (e.g., reduced bending stiffness)
- record `method=AISC360_DAM` and `dam_ok=True/False` in results

Deliverable: a single method call like:

- `LoadedFrame(..., settings=FrameAnalysisSettings(method="AISC360_DAM"))`

### 4.3 EC3 wrapper

Wrapper responsibilities:

- decide when second-order is required (via stability metric / user override)
- ensure global imperfections are included (notional loads or initial sway)
- apply EC3 stiffness conventions (kept configurable; don’t hardcode too early)

Deliverable: `FrameAnalysisSettings(method="EC3_GLOBAL_IMPERFECTIONS")` produces a reproducible analysis recipe.

### 4.4 AS 4100 wrapper

Wrapper responsibilities:

- decide amplified-first-order vs second-order (threshold logic; user override)
- implement the AS notional horizontal force convention via `notional_factor`

Deliverable: `FrameAnalysisSettings(method="AS4100_SECOND_ORDER")`.

---

## 5) Results, envelopes, and checker integration

### 5.1 Results tagging

Persist to `LoadedFrame` results:

- `analysis_method`
- `second_order_enabled`
- `imperfection_model`
- `stiffness_rules`
- `converged`, `iterations`

Add two checker-safety fields:

- `design_grade_ok: bool` (stricter than `converged`)
- `design_grade_notes: list[str]` (reasons if false: nonconvergence, stabilization DOFs added, stepping hit max, etc.)

### 5.2 Member demand extraction should not change

Keep the existing “recover member end forces from global solution” approach.
Second-order only changes the equilibrium solution, not the extraction API.

#### 5.2.1 Define a canonical member demand API

Standardize a single “member demand” payload used by *all* checkers:

- member end forces in **local axes**: $(N, V_y, V_z, T, M_y, M_z)$ at both ends
- internal action envelopes (or station samples) along member for the same set
- explicit sign conventions documented once and referenced everywhere

This prevents checkers from re-interpreting signs or recomputing demands inconsistently.

### 5.3 AISC checks should enforce consistency

Example policy for Chapter E:

- If `analysis_method == AISC360_DAM` and converged: default `K=1.0` (unless user overrides)
- Otherwise: require explicit `K` or provide a non-DAM path

EC3/AS checkers can later use the same results metadata.

#### 5.3.1 Numerical stabilization must never imply physical restraint

Auto-added constraints (e.g., fixing RX/RY/RZ at axial-only nodes) must be:

- reported distinctly as **numerical stabilization**,
- excluded from any “physical restraint / bracing” inference used by design checks.

---

## 6) Verification plan (small deterministic benchmarks)

Add a few “known behavior” checks (not necessarily exact closed-form):

1. **Sway frame under gravity + lateral**: second-order should increase drift and moments.
2. **Leaned column (P–Δ)**: compare to classic amplification trend (qualitative + approximate).
3. **Symmetry**: symmetric frame + symmetric loads stays symmetric in results.

Also re-run existing scripts (e.g., tf/tf.py) and compare:

- linear vs second-order drift
- member utilization changes (Chapter E especially)

---

## 7) Delivery phases (suggested sequencing)

1. **Plumbing only**: `FrameAnalysisSettings`, `analysis_method` tagging, no behavior change.
2. **Notional loads**: generate + apply, still linear.
3. **Stiffness modifiers**: apply reduced EI via settings.
4. **Second-order loop**: geometric stiffness + iteration + convergence.
5. **Wrappers**: AISC DAM preset, EC3 preset, AS preset.
6. **Node springs (6-DOF)**: improve restraint realism across all codes.

7. **Future-proof springs**: design API to eventually accept a full 6×6 coupled spring matrix per node (even if initially only diagonal terms are implemented).

