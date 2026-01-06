A workable structure is: keep **model objects**, **loading**, **analysis/solver**, **results**, **plotting**, and **design checks** in separate modules. Don’t put math/solver code in `member.py`/`frame.py`; those should stay as data + validation + convenience wrappers.

Here’s a concrete package layout and what belongs where.

---

## Proposed package structure

```
beamy/
  pyproject.toml
  README.md
  src/beamy/
    __init__.py

    model/
      __init__.py
      node.py
      member.py
      frame.py
      dof.py
      connections.py

    loads/
      __init__.py
      loadcase.py
      nodal_loads.py
      member_loads.py
      supports.py

    analysis/
      __init__.py
      model.py
      settings.py
      mesh.py
      assembly.py
      solver_linear.py
      solver_second_order.py
      recovery.py
      transformations.py
      elements/
        __init__.py
        beam3d.py
        truss3d.py
        cable3d.py


    results/
      __init__.py
      frame_result.py
      member_result.py

    viz/
      __init__.py
      plot_frame.py
      plot_member.py

    design/
      __init__.py
      aisc/
        __init__.py
        compression.py
        flexure.py
        interaction.py

    io/
      __init__.py
      json.py
      csv.py              # optional exports

```

* **pyproject.toml**: packaging metadata, dependencies (numpy, scipy optional, matplotlib), dev deps (pytest).
* **README.md**: minimal examples (build frame, apply loadcase, analyze, plot, design check).

---

```
src/beamy/
  __init__.py
```

* Re-export the public API: `Frame`, `Member`, `Node`, `LoadCase`, `AnalysisSettings`.
* Keep it thin (no logic).

---

## model/  (user-facing structural model; no solver state)

```
model/
  __init__.py
```

* Re-export key model classes: `Node`, `Member`, `Frame`, connection refs.

```
model/node.py
```

* `Node(id, xyz, meta)`
* Basic validation (`xyz` shape, numeric types).
* No tolerance/merging logic (that belongs in `Frame`).

```
model/member.py
```

* `Member` dataclass: endpoints (`end_i`, `end_j`), `kind`, properties refs, releases, stations, second-order flags.
* **DesignProperties**: Separate physical length from design parameters.
  * `Lx`, `Ly` (unbraced lengths)
  * `Kx`, `Ky`, `Kz` (effective length factors)
  * `Cb` (moment gradient factor)
* Derived geometry helpers that *require* a `Frame` (e.g., `length(frame)`, `local_axes(frame)`).
* No global assembly, no solving, no caching results.

```
model/frame.py
```

* `Frame` dataclass holding `nodes`, `members`, `connections`, `merge_tol`.
* **Registries**: `materials: dict[str, Material]`, `sections: SectionRegistry`. 
  * **SectionRegistry**: Smart lookup that auto-populates standard shapes (e.g., "W14X90") from a DB/CSV on first access.
  * Ensures single source of truth and validation at build time.
  * Members store `material_id` / `section_id` strings, resolved against these registries.
* Model-building methods: `add_node`, `add_member`, `connect_end_to_member`, optional `remove_*`.
* `validate()` for IDs, missing properties, invalid station ranges, etc.
* Convenience wrappers: `analyze(loadcase, settings)` and `plot(...)`.
  * **Import Rule**: Use local imports (e.g., `from beamy.analysis.engine import run_analysis`) inside these methods to avoid circular dependencies.

```
model/dof.py
```

* DOF conventions and masks.

  * Define DOF order for 3D frame: `[UX, UY, UZ, RX, RY, RZ]` (and optionally 2D later).
  * Define two mask types: `ReleaseMask6` (True=released) vs `RestraintMask6` (True=restrained).
* Helpers: `dof_index("RY")`, mask constructors, validation utilities.

```
model/connections.py
```

* Typed connection objects:

  * `EndRef(member_id, end="i|j")`
  * `StationRef(member_id, s)`
  * `EndToStation(end_ref, station_ref, dofs=None)`
  * (later) MPC / rigid link definitions
* Validation helpers (e.g., `0 <= s <= 1`).
* No meshing/splitting logic (that’s `analysis/mesh.py`).

---

## loads/  (inputs per analysis case; no solver math)

```
loads/
  __init__.py
```

* Re-export `LoadCase`, load/support classes.

```
loads/loadcase.py
```

* `LoadCase(name, nodal_loads, member_loads, supports, prescribed_displacements, springs)`
* Lightweight validation (referenced ids exist if a Frame is provided to validate against).
* No conversion into global force vectors (that’s `analysis/assembly.py`).

```
loads/nodal_loads.py
```

* Data types for nodal forces/moments:

  * `NodalLoad(node_id, forces=(Fx,Fy,Fz), moments=(Mx,My,Mz), coord_sys="global|local")`
* Helpers for combining loads or converting to a standard representation.

```
loads/member_loads.py
```

* Member-applied loads:

  * `MemberPointLoad(member_id, s, vector, coord_sys)`
  * `MemberPointMoment(member_id, s, vector, coord_sys)`
  * `MemberDistributedLoad(member_id, s0, s1, w0, w1, coord_sys)`
* Only stores definitions; no equivalent nodal load calculations here (put that in `analysis/recovery.py` or `analysis/assembly.py` depending on approach).

```
loads/supports.py
```

* Boundary conditions and constraints-to-ground:

  * `Support(target=NodeId|StationRef, restrained_dofs=RestraintMask6)`
  * (optional) `Spring(target, k_trans, k_rot)`
  * `PrescribedDisplacement(target, values)`
* Again: just definitions.

---

## analysis/  (all FE mechanics, meshing, assembly, solving, recovery)

```
analysis/
  __init__.py
```

* Export the main entry point, e.g. `run_analysis(frame, loadcase, settings)`.

```
analysis/model.py
```

* Internal solver model structures (created inside analyze):

  * `AnalysisNode(id, xyz)`
  * `Element(id, type, node_i, node_j, parent_member_id, local_axes, releases_end_flags, etc.)`
  * Mappings: `member_to_elements`, `station_node_map`
  * Per-element state containers for second-order (axial force history, corotational frame, etc.)
* This file should be “dumb containers + typing”, not algorithms.

```
analysis/settings.py
```

* `AnalysisSettings`:

  * `dimension: "3d"` (or future `"2d"`)
  * `solver: "linear" | "second_order"`
  * `formulation: "linear" | "pdelta" | "corotational"`
  * tolerances, max iters, load steps, convergence criteria
  * options for automatic station insertion / mesh density

```
analysis/mesh.py
```

* Builds `AnalysisModel` from `Frame + LoadCase`:
  * **Automated Station Discovery (Station Manager)**:
    * **Collect** all $s$ values from `Member.stations`, `PointLoad`, `Support`, and `Connection`.
    * **Deduplicate** using `math.isclose` with tolerance (e.g., `merge_tol/L`).
    * **Sort** `[0.0, ..., 1.0]`.
    * **P-Delta Refinement**: If `formulation="pdelta"`, force insertion of at least one mid-point station ($s=0.5$) to capture member curvature.
  * create analysis nodes at stations
  * split members into analysis elements
  * apply end-to-station connectivity (node sharing) by wiring endpoints to the inserted station node
* Output: `AnalysisModel`

```
analysis/transformations.py
```

* Local/global coordinate transforms:
  * direction cosines
  * rotation matrices
  * building element transformation `T`
  * **Vertical Member Handling**: 
    * Detects if member is parallel to global up vector.
    * Falls back to "Global North" or "Global East" reference to resolve orientation ambiguity.
  * (second-order) updated orientation/corotational helpers if you implement them

```
analysis/elements/__init__.py
```

* Element factory/registry:

  * `make_element(element, frame, settings)` or dispatch table by type.

```
analysis/elements/beam3d.py
```

* Beam element routines (core structural element):
  * local stiffness `k_material(E, G, A, Iy, Iz, J, L, ...)`
  * geometric stiffness `k_geometric(N, L, ...)` (for p-delta)
  * fixed-end/equivalent nodal loads for member loads (if you implement here)
  * **Static Condensation**:
    * Handle releases at the **Element Level**.
    * Partition $K_e$ and condense out released DOFs before returning to assembler.
    * Store **Condensation Matrix** ($K_{cc}^{-1} K_{cr}$) for instant internal DOF recovery.
  * recovery: end forces from `u_local`, section forces along stations
* Keep functions testable and mostly pure.

```
analysis/elements/truss3d.py
```

* Truss element stiffness/recovery:

  * axial-only stiffness
  * axial force recovery
* Optionally geometric stiffness for second-order if you support truss P-Δ.

```
analysis/elements/cable3d.py
```

* Cable/truss-with-tension-only behavior (if you implement it):

  * tension-only iteration logic likely belongs in solver (or element “active/inactive” flag)
  * local stiffness same as truss when active
* Keep it minimal early; cables add nonlinearity.

```
analysis/assembly.py
```

* Global system building:
  * DOF numbering (map each analysis node DOF to equation index)
  * assemble global `K` and `F`
  * apply releases: Already handled by `Beam3D` returning condensed matrices.
  * apply supports/constraints (restraints, prescribed displacements, springs)
* Output: assembled system (and bookkeeping to compute reactions)

```
analysis/solver_linear.py
```

* Solve one linear system:

  * `u = solve(K, F)`
  * compute reactions
* Should not care about plotting/results formatting.

```
analysis/solver_second_order.py
```

* Second-order analysis driver:

  * load stepping + Newton/modified Newton
  * per-iteration: assemble tangent stiffness `K_t = K_m + K_g`
  * update element axial forces and states used for `K_g`
  * convergence checks, line search (optional)
* Produces final `u` plus iteration history.

```
analysis/recovery.py
```

* Post-processing into engineering outputs:

  * convert `u` into element end forces
  * compute internal force diagrams at requested stations
  * compute stresses (axial + bending) per station
  * aggregate element results back to parent members
* Produces `FrameResult`/`MemberResult` structures.

---

## results/  (structured outputs + query helpers)

```
results/
  __init__.py
```

* Re-export `FrameResult`, `MemberResult`.

```
results/frame_result.py
```

* `FrameResult` dataclass:
  * **Standalone Snapshot**: Indexed by element IDs at time of solve. Does not reference mutable Frame state.
  * node displacements (by original node id and/or analysis node id)
  * reactions at supports
  * per-element end forces
  * per-member diagrams + stress results
  * solver metadata: convergence, steps/iters
* Keep it data-first; include small convenience accessors.

```
results/member_result.py
```

* `MemberResult` dataclass:

  * diagrams (N/V/M/T), deflection shape, stress arrays
  * mapping from station `s` → values
* Typically created from `FrameResult.member(member_id)`.

(If you later add `results/extract.py`, it can hold table exports/envelopes.)

---

## viz/  (matplotlib-only; read-only)

```
viz/
  __init__.py
```

* Export plotting functions.

```
viz/plot_frame.py
```

* Plot undeformed frame:

  * nodes, members, labels, supports markers
  * optional: show local axes
* Input: `Frame` (and maybe style options). No solver math.

```
viz/plot_member.py
```

* Member-specific plots:

  * diagrams (N/V/M/T)
  * stress plots
  * deflection along member
* Input: `Frame`, `MemberResult` or `FrameResult`.

(You may later prefer `model.py/deformed.py/diagrams.py`, but this is fine for now.)

---

## design/  (code checks; independent from solver internals)

```
design/
  __init__.py
```

* Re-export main design entry points.

```
design/aisc/__init__.py
```

* Public API for AISC checks, e.g. `check_member_aisc360(...)`.

```
design/aisc/compression.py
```

* Chapter E functions (compression/buckling):

  * compute `Pn`, `phiPn`/`OmegaPn` depending on method
  * requires section/material props + effective lengths

```
design/aisc/flexure.py
```

* Chapter F functions (bending strength):

  * compute `Mn` for strong/weak axis
  * LTB checks need `Lb`, `Cb`, etc.

```
design/aisc/interaction.py
```

* Chapter H (interaction equations):

  * combine `Pu`, `Mux`, `Muy` with capacities from E/F
  * returns utilization ratio + governing equation

Design should accept:

* demands from `MemberResult` / `FrameResult`
* plus user-provided design parameters (K, Lb, Cb, bracing, etc.)

---

## io/  (import/export)

```
io/
  __init__.py
```

* Re-export `load_json`, `save_json`, etc.

```
io/json.py
```

* Serialize/deserialize:

  * `Frame`, `LoadCase`, `AnalysisSettings`
  * careful with dataclasses and versioning
* Keep stable schema (include a version field).

```
io/csv.py
```

* Export result tables:

  * nodal displacements
  * reactions
  * member force/stress summaries
* No plotting, no solving.

---

### One naming note (to avoid bugs)

In `model/dof.py`, keep these separate:

* **ReleaseMask** (member end; True=released)
* **RestraintMask** (supports; True=restrained)

Same shape, different semantics.

---

If you want, I can also propose the *imports direction* rule (who is allowed to import whom) to prevent circular imports—this matters a lot once `frame.py` has convenience methods calling `analysis`.

---

## Architecture Rules

### 1. The "Import Direction" Rule

To prevent circular dependencies, adhere to this one-way dependency stack:

1.  **Level 0 (Base)**: `dof.py`, `nodes.py`, `materials.py` (no imports from other levels).
2.  **Level 1 (Entities)**: `member.py`, `connections.py` (imports Level 0).
3.  **Level 2 (Container)**: `frame.py` (imports Level 0 & 1). **Never** imports `analysis/` or `results/` at module level.
4.  **Level 3 (Scenario)**: `loadcase.py`, `supports.py` (imports Level 2).
5.  **Level 4 (Mechanics)**: `analysis/` (imports Frame, LoadCase).
6.  **Level 5 (Output)**: `results/`, `viz/`, `design/` (imports all below).

**Injection/Delegation**: If `Frame` needs to `analyze()`, use a local import inside the method:
```python
def analyze(self, ...):
    from beamy.analysis import run_analysis
    return run_analysis(self, ...)
```
