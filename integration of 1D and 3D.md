
# Integration of 1D and 3D (Redesign Plan)

This note defines how **1D member analysis** and **3D frame analysis** should converge on a single, consistent internal workflow.

## Motivation

The frame solver already has a robust mechanism that takes a **loaded frame** and produces a **continuous internal-action profile** for a particular member:

- `LoadedFrame.demand_provider.actions(member_id, points=...)`
- `LoadedFrame.demand_provider.actions_original(original_member_id, ...)`

This path uses **global frame analysis end forces** and reconstructs actions via **equilibrium** (`MemberResultsDirect` / `_compute_internal_forces_direct`) rather than re-solving a member in isolation with artificial boundary conditions.

Goal: make “1D analysis” consume the same pipeline so that 1D and 3D results (and design checks) are consistent by construction.

## Naming decision

### `LoadedMember` as the canonical name

The object representing "beam + loads + analysis results" is now `LoadedMember` throughout the codebase.

Rationale: many "members" are not strictly beams (e.g., frame elements, segments after splitting), and the future architecture uses the same extraction concept whether analyzing a standalone member or extracting from a frame.

## Chosen architecture: Strategy A (recommended)

### Strategy A — 1D delegates to frame analysis internally

The 1D path should build and solve a minimal 3D model under the hood:

1. Build a single-member `Frame` representing the 1D member.
2. Insert nodes at all relevant x-stations:
	 - support locations
	 - point load locations
	 - moment locations
	 - distributed load start/end locations
3. Convert `LoadCase` → `FrameLoadCase` by attaching loads to that single member:
	 - `PointForce` → `MemberPointForce`
	 - `Moment` → `MemberPointMoment`
	 - `DistributedForce` → `MemberDistributedForce`
4. Apply supports by converting each `Support(x, type)` into node fixities on the appropriate inserted node.
5. Run `LoadedFrame.analyze()`.
6. Query the member’s action profile via `demand_provider.actions_original(...)`.

Outputs:

- A canonical `MemberActionProfile` containing continuous arrays for $(N, V_y, V_z, T, M_y, M_z)$ along $x$.
- Optional: displacements can either be extracted from the single-member frame solution or remain a 1D-only feature until unified.

Why Strategy A:

- One solver for end forces and reactions → less divergence over time.
- Consistent sign conventions and stiffness formulation.
- Naturally supports the same concepts as frames (releases, springs, second-order options later).

## Shared member-demand engine (critical refactor)

Regardless of Strategy A vs B, we want a single “member line demand” implementation.

Plan:

- Move equilibrium reconstruction into a shared module (e.g., `src/beamy/core/member_line.py`).
	- Input: member length, local end forces, distributed line loads (local coords)
	- Output: `MemberActionProfile`
- Make `MemberDemandProvider` call this shared implementation.
- Make the 1D wrapper (`LoadedMember`) call the same implementation.

This guarantees that plots/checks see identical demand profiles for a given set of end forces + member loads.

## API changes (proposed)

### New/primary

- `LoadedMember(member: Beam1D, loads: LoadCase, *, backend='frame', settings=None)`
	- `backend='frame'` uses Strategy A
	- `backend='legacy_1d'` temporarily keeps the old 1D solver for comparison/testing

Key properties/methods:

- `.actions(points=...) -> MemberActionProfile`
- `.envelopes(points=...) -> dict`
- `.plot_*` and AISC checks consume `MemberActionProfile`

## Validation plan

1. **Equivalence tests**: single-member frame vs `LoadedMember(backend='frame')` (actions + reactions).
2. **Classical benchmarks**: cantilever end load, simply supported point load, UDL.
3. **Regression**: ensure frame member extraction (via `demand_provider`) matches prior verified outputs.

## Notes / open questions

- Displacements: do we want `LoadedMember.deflection()` to remain supported immediately? If so, it should pull from the frame nodal displacements for the single-member model and interpolate similarly to the frame plotting shape functions.
- Coordinate conventions: define one documented convention for local axes and action signs and enforce it everywhere (1D + 3D + checks).

