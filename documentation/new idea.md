
## New API + architecture idea (frame-first, member-as-frame)

This note fleshes out the “everything is a frame” strategy for Beamy.


### Key principle: members *can* have loads, but they’re rebuilt on solve

Members should be allowed to carry loads (it’s a very natural user workflow), **but it must be clear what that means**:

- **`LoadCase` remains the canonical input container** for a frame solve.
- **`Member.loads` is treated as a derived / cached “applied-to-this-member” view** that the solver can rebuild.

So the rule is explicit:

- **When you solve a frame, you update the loads applied to each member.**
  - Any prior “applied loads on members” are cleared/overwritten.
  - Loads are re-derived from the active `LoadCase` (including any solver rewriting due to auto-splitting).

This keeps the API ergonomic (“members have loads”) while keeping the analysis deterministic (“solve defines the applied loads used”).

---

### Core objects (conceptual model)

#### `Member` (geometry + element behavior)

A `Member` is one finite element "line" with:

- **Identity**: `id`
- **Geometry**: `start` and `end` positions (3D coordinate arrays)
- **Section/material/orientation**
- **Element type**: `"beam" | "truss" | "cable"`
- **End releases**: internal connection flexibility (hinges) at member ends
- **Applied loads**: none by default (List of PointForces, PointMoments, DistributedForce objects)

The following is properties of the Member class that are generated when a frame is analysed
- **Loads**: none by default (List of PointForces, PointMoments, DistributedForce objects)
This is generated when the frame is analysed
It would obviously include the applied loads, and loads where it connects to other members

- **End constraints** (optional): *additional* fixed DOFs at member ends (a modeling convenience)
- **Point supports** (optional): supports along the member geometry (restraints-to-ground at interior points)
These 2 are basically the same and should be part of the frame


Important distinction:

- **Releases** change the *member stiffness* (joint/connection modeling).
- **Constraints** change the *global boundary conditions* (restraint-to-ground at member ends).






#### `Frame` (the structure)

`Frame` is just:

- `members: list[Member]` (input: defined by start/end positions)
- `nodes: dict[Node]` (auto-generated from member endpoints) (not exactly sure what is required here)

No loads. No results. It's the *model*.

**Construction**:
- `Frame(members)` - auto-generates nodes from member endpoints


#### `LoadCase` (applied loads and boundary conditions container)

LoadCase is made of loads and supports
Loads and supports may be applied to nodes or members

`LoadCase` stores **applied** actions and **boundary conditions**:

- Nodal loads: forces/moments/springs
- Member loads: point forces/moments, distributed loads
- Member point supports: restraints-to-ground at specific positions along members
- Member supports: restraints-to-ground applied to all nodes along a member

After analysis you ask for **demand** (internal actions, deflections, stresses).

Important distinction:

- **Applied loads**: what was placed on nodes/members (from `LoadCase`, then distributed/rewritten during solve)
- **Demand**: the solved internal response, which can exist on any member **even if it has no direct applied loads**
  due to global compatibility/load-path.
- **Point supports**: boundary conditions (part of `LoadCase`), NOT part of member geometry.

#### `Frame.analyze(...)` (solved state for one load case)

`frame.analyze(load_case, settings)`:

- Validates the inputs (unknown node/member ids, invalid local axes references, etc.)
- Expands the model for attachment points (see below)
- Assembles/solves/recover results
- Caches a `FrameAnalysisResult` on the `Frame`
- Exposes `frame.member_demand(member_id)` and plotting helpers (viz functions accept `Frame`)

#### `LoadedMember` (a convenience wrapper)

`LoadedMember` is “beam analysis” but implemented as:

> Build a 2-node `Frame` + one `Member`, apply end supports, then run `frame.analyze(load_case)`.

This guarantees that “beam” and “frame” share the same solver and result interpretation.

---

### The tricky but important part: attachment points & auto-splitting

Users want to apply:

- A point load at x = 2.3 m
- A spring support at midspan
- A distributed load that starts/stops inside the member
- Two members that intersect at a node along a “continuous” member

These should “just work” without forcing users to manually split members.

**Strategy** (solver-internal):

- **Insert real nodes** at any “attachment point” along a member:
  - member point forces / point moments
  - distributed load endpoints
  - member point supports (restraints-to-ground along the member geometry)
  - other frame nodes that lie on the member line (connectivity)
- **Split the member** into solver segments between those nodes.
- **Rewrite loads** so the solver only sees:
  - nodal forces/moments/springs
  - distributed loads that live fully on a solver segment
- Preserve a mapping “original member → ordered solver segments” so demand can be queried
  as if the original member was continuous.

This is already implemented in the solver expansion step in `src/beamy/frame/analysis.py`:
it expands the model, solves it, then provides “bundled” member demand back to the user via `frame.member_demand(...)`.

---

### Boundary conditions: supports vs releases vs constraints (practical semantics)

Beamy uses 6 DOFs per node:

- `[UX, UY, UZ, RX, RY, RZ]`

#### Node support (`Node.support`: 6 digits)

Represents restraint-to-ground at a node:

- `"1"` means fixed, `"0"` means free.
- Examples:
  - `"111111"` fixed support
  - `"111000"` pin (translations fixed, rotations free)

#### Member end releases (`Member.releases`: 12 digits)

Represents connection flexibility at member ends (hinges/releases), by modifying element stiffness:

- 12 digits: first 6 = start, last 6 = end (same DOF ordering).
- `"1"` means released (free), `"0"` means connected (rigid).
- Releases do **not** represent restraint-to-ground.

#### Member end constraints (`Member.constraints`: 12 digits, optional)

Represents “treat this member end DOF as fixed to ground” (a modeling convenience).

If we keep this concept, the rule should be:

- Constraints behave like additional support bits applied at the member’s start/end node.
- They are additive with node supports.

#### Member point supports (`Member.point_supports`: supports along geometry)

Represents physical restraints-to-ground at an interior point on the member without the user
manually splitting the member.

Internally this requires inserting a node at that x-location and applying a node support there.

---

### Solve pipeline (what happens when you “analyze”)

#### Input

- A `Frame` (nodes + members)
- A `LoadCase` (applied nodal/member loads)
- Optional `FrameAnalysisSettings`

#### Steps

- Validate references (node ids, member ids, local-axis references).
- Expand model for attachment points (auto-split, load rewriting).
- Assemble global stiffness, apply releases, apply supports/constraints.
- Solve global displacements.
- Recover member end forces (local).
- Provide:
  - nodal displacements
  - reactions (with careful separation of physical vs numerical stabilization reactions)
  - member demand access (actions/deflections/stresses sampled along the member)

#### Output

One analysis produces one “solved state” container (`FrameAnalysisResult`) plus demand accessors.
The important ergonomic point: **results are retrieved by queries**, not by mutating the model.

---

### Plotting strategy

Plotting should hang off analysis objects (or accept them), because plots are driven by results:

- `plot_frame(frame, load_case=None, save_path="...svg")`
- `plot_deflection(analysis, save_path="...svg")`
- `plot_von_mises(analysis, save_path="...svg")`
- `plot_member_diagrams(analysis, member_id="M1", save_path="...svg")`

Key rule: **plot saving should default to SVG** (and examples/docs should show `.svg`).

---

### “Analyze a member in isolation” (member → frame)

This is the user-facing convenience that makes the API feel “beam-first” while staying
frame-consistent:

- Create two nodes at the member ends.
- Apply end supports on those nodes.
- Create a `Frame` with one member.
- Run `frame.analyze(load_case)`.

The only special-case behavior should be validation (e.g., a standalone member should not accept
loads that target unknown nodes or other members).

---

### Where this simplifies the codebase

If we enforce these boundaries:

- Geometry is immutable-ish (or at least separate from loads/results).
- `LoadCase` is the only applied load container.
- The solver owns splitting/rewriting logic.
- Demand/recovery lives in one place and is consistent for beams + frames.

Then:

- We can delete adapters that “re-solve members as 1D beams” to approximate internal forces.
- Checks (AISC, etc.) can run directly on **solved member demand** (actions profiles), not on raw loads.
- Plotting becomes consistent across “beam” and “frame”.

---

### Suggested public API shape (high-level)

Keep the public surface small and composable:

- **Model**
  - `Node`, `Member`, `Frame`
- **Loads**
  - `LoadCase` + load dataclasses (nodal + member)
- **Analysis**
  - `frame.analyze(load_case, settings)` (mutates the *Frame's cached results*, not the geometry)
  - `LoadedMember(..., load_case, support_start, support_end, ...)` (wrapper)
- **Results**
  - `FrameAnalysisResult`
  - `MemberDemand` / `MemberActionProfile` (query objects)
- **Viz**
  - pure functions that accept `Frame` (optionally `load_case` for geometry plots) and save `.svg`

---

### Migration notes (if refactoring from an older member-load approach)

- Stop mutating members when solving (no “clear member loads”).
- Move all applied loading into `LoadCase`.
- Provide a clear “analyze once per load case” call (`frame.analyze(load_case)`).
- Keep backwards-compatible helpers temporarily (thin wrappers that build `LoadCase` + call analysis).
- Add deprecation warnings for any API that stores loads on members.

