

# Core model (user-facing)

## Node

Keep Nodes simple, but add a tolerance-friendly identity strategy at Frame level.

**Attributes**

* `id: NodeId`
* `xyz: Vec3`
* `meta: dict = {}` (optional)

---

## Member (physical member, user-facing)

### A) Identity + endpoints

* `id: MemberId`
* `end_i: NodeId`
* `end_j: NodeId`
* `kind: Literal["beam","truss","cable"]`
* `group: str | None = None`
* `meta: dict = {}`


### B) Properties (structural)

* `material: MaterialRef` (resolved from Frame registry)
* `section: SectionRef` (resolved from Frame registry)
* `orientation: OrientationDef`
  One of:

  * `local_z: Vec3`
  * `roll_angle: float`
  * `ref_node: NodeId`
* **DesignProperties**: Separate physical length from design parameters.
  * `Lx`, `Ly` (unbraced lengths)
  * `Kx`, `Ky`, `Kz` (effective length factors)
  * `Cb` (moment gradient factor)

### C) End behavior (connection stiffness modifiers)

* `releases_i: DofReleaseMask`
* `releases_j: DofReleaseMask`

**Optional but recommended** (future-proof)

* `end_offsets_i: Vec3 = (0,0,0)`
* `end_offsets_j: Vec3 = (0,0,0)`
  (Offsets help later if you implement rigid links, eccentricity, panel zones.)

### D) Meshing / stations (not loads)

Replace generic `attachments` with a clearer concept: **stations**.

* `stations: set[float] = {}`  (normalized 0..1)

Use stations for:

* “must split here if anything happens here”
* user mesh seeds
* connection landing points (Frame will add these when you connect end-to-interior)

This avoids mixing “topology facts” with “mesh constraints”. A station is just: “this member may need a node here”.

### E) Analysis behavior (second-order readiness)

Keep it minimal but explicit:

* `formulation: Literal["linear","pdelta","corotational"] = "pdelta"`
* `integration_points: int = 2` (optional)

### F) Convenience (member-only analysis)

Keep these orchestration helpers, but don’t store state:

* `to_subframe(frame: Frame | None = None) -> Frame`
* `analyze(load_case: LoadCase, settings: AnalysisSettings) -> MemberResult`

---

# Frame (structure/topology, user-facing)

## Frame attributes

### A) Topology

* `members: dict[MemberId, Member]`
* `nodes: dict[NodeId, Node]`
* **Registries**: `materials: dict[str, Material]`, `sections: dict[str, Section]`

**Recommended workflow (explicit)**

* Define materials and sections in the registry.
* Create nodes first (`Frame.add_node(...)`), then create members referencing those node ids, then add members to the frame.
* After all members (and any explicit connections) are added, call a `Frame.build()` method (or do this automatically at the start of `analyze()`).

**What `Frame.build()` should do**

* Validate references (each `Member.end_i/end_j` exists in `Frame.nodes`, no duplicate ids, material/section IDs exist in registry, etc.).
* Create implicit connectivity: member ends that share the same `NodeId` are a joint (no extra “connection object” required).
* Materialize split requirements via **Automated Station Discovery**:
  * for every `EndToStation(... target=StationRef(member_id, s) ...)`, ensure `s` exists in that member’s `Member.stations`
  * for every load/support targeting a `StationRef`, ensure the station exists too
  * Point load locations.
* (Optionally) normalize/clean stations (dedupe, clamp to 0..1, drop ~0/~1 if you don’t want end-stations).

### B) Connectivity constraints

Instead of a generic `connections: list[Connection]` with a `data: dict`, define a typed connection representation:

**EndRef**

* `(member_id, end: Literal["i","j"])`

**StationRef**

* `(member_id, s: float)`  (0..1)

**Connection types**

1. `NodeJoint(node: NodeId)` — implicit when multiple member ends share same node id
2. `EndToStation(end: EndRef, target: StationRef, dofs: DofMask)` — end connects into middle of other member
3. (optional) `MPC(...)` for rigid links, equal DOFs, etc.

So Frame stores:

* `connections: list[EndToStation | MPC | ...]`

Why this helps:

* You can model “member end connects to middle” *without* inventing fake nodes prematurely.
* `analyze()` will insert the needed node at that station and split the host member.

### C) Node generation policy

* `merge_tol: float = 1e-6`
* `node_id_strategy: Literal["explicit","auto_merge"] = "auto_merge"`
* (optional) spatial hashing config if you want performant merging

### D) Public methods

* `add_node(xyz, id=None) -> NodeId`
* `add_member(end_i, end_j, ...) -> MemberId`
* `connect_end_to_member(end: EndRef, target: StationRef, dofs="all")`
* `validate()`
* `plot(model=True, result=None, deformed=False, ...)`
* `analyze(load_case, settings) -> FrameResult`

---

# Applied vs resultant (unchanged, but tighten definitions)

## LoadCase (input only)

Don’t put loads/supports on Member/Frame.

**Recommended schema**

* `nodal_loads: list[NodalLoad(node_id, dof_vector)]`
* `member_point_loads: list[MemberPointLoad(member_id, s, local_or_global, vector)]`
* `member_dist_loads: list[MemberDistLoad(member_id, s0, s1, ...) ]`
* `supports: list[Support(target: NodeId | StationRef, restrained_dofs)]`
* `prescribed_displacements: list[...]`
* `springs: list[...]`

Note: supports can target a StationRef; analysis will create the node and apply the support there.

## Results (output only)

* `FrameResult` includes:
  * **Standalone Snapshot**: Indexed by element IDs at time of solve.
  * `u: dict[AnalysisNodeId, dof_values]`
  * `reactions: dict[AnalysisNodeId, dof_values]`
  * `element_end_forces: dict[ElementId, EndForces]`
  * `member_diagrams: dict[MemberId, DiagramData]` (aggregated from element pieces)
  * nonlinear: `converged`, `step_history`, `iter_history`

---

# Solver-facing (internal, created in analyze)

This is where second-order needs extra structure.

## AnalysisModel (internal)

* `analysis_nodes: dict[ANodeId, Node]`
* `elements: dict[ElementId, Element]`
* `member_to_elements: dict[MemberId, list[ElementId]]`
* `station_node_map: dict[(MemberId, s), ANodeId]`
* `element_state: dict[ElementId, ElementState]` (for corotational / second-order)

## Element (internal)

An Element is a *piece* of a Member after splitting at stations:

* `id: ElementId`
* `parent_member: MemberId`
* `node_i: ANodeId`
* `node_j: ANodeId`
* `releases_i / releases_j` (inherited if this piece touches the original end)
* `local_axes` / orientation
* methods:

  * `k_material(...)`
  * `k_geometric(N, ...)`
  * `recover_forces(u, ...)`

This keeps second-order out of the user model and makes your nonlinear implementation sane.

---

# What changes vs your current version (delta)

### 1) Replace `attachments` with two clearer concepts

* **Topology connections** live in `Frame.connections` as typed objects (EndToStation, etc.)
* **Meshing stations** live in `Member.stations` (a set of `s` positions)

During `Frame.connect_end_to_member(...)`, you record the intent:

* add a `Connection(EndToStation(...))`

Then `Frame.build()` / `analyze()` materializes the split requirement:

* ensure the target member has station `s` in `Member.stations` (so analysis knows it must split there)

### 2) Remove `Connection.data: dict`

Use typed fields. You’ll thank yourself later when validating and debugging.

### 3) Add second-order readiness to Member

Minimal flag: `formulation` (e.g. `"linear"`, `"pdelta"`, `"corotational"`).

### 4) Keep Frame/Member state immutable-ish

No solved forces stored on Member. Everything solved goes into results.

---

# Improved Python skeleton (tight, typed)

```python
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Tuple, Set, Union

Vec3 = Tuple[float, float, float]
NodeId = str
MemberId = str

# Convention (3D frame):
# DOF order = [UX, UY, UZ, RX, RY, RZ]
Mask6 = Tuple[bool, bool, bool, bool, bool, bool]
ReleaseMask6 = Mask6     # member end releases: True = released (no stiffness transmitted)
RestraintMask6 = Mask6   # supports/restraints: True = restrained (fixed to ground)

# --- Core ---

@dataclass(frozen=True)
class Node:
    id: NodeId
    xyz: Vec3
    meta: dict = field(default_factory=dict)

@dataclass(frozen=True)
class Orientation:
    local_z: Optional[Vec3] = None
    roll_angle: Optional[float] = None
    ref_node: Optional[NodeId] = None
    # If vertical, solver will fall back to Global North/East

@dataclass(frozen=True)
class DesignProperties:
    Lx: Optional[float] = None
    Ly: Optional[float] = None
    Kx: float = 1.0
    Ky: float = 1.0
    Cb: float = 1.0

@dataclass(frozen=True)
class Member:
    id: MemberId
    end_i: NodeId
    end_j: NodeId
    kind: Literal["beam", "truss", "cable"]
    material_id: str
    section_id: str
    orientation: Orientation = Orientation()
    design: DesignProperties = DesignProperties()
    releases_i: ReleaseMask6 = (False, False, False, False, False, False)
    releases_j: ReleaseMask6 = (False, False, False, False, False, False)

    # meshing stations (0..1), no loads here.
    # Store as a unique, sorted tuple for immutability + stable behavior.
    stations: Tuple[float, ...] = ()

    # second-order controls
    formulation: Literal["linear", "pdelta", "corotational"] = "pdelta"
    integration_points: int = 2

    group: Optional[str] = None
    meta: dict = field(default_factory=dict)

@dataclass(frozen=True)
class EndRef:
    member_id: MemberId
    end: Literal["i", "j"]

@dataclass(frozen=True)
class StationRef:
    member_id: MemberId
    s: float  # 0..1

@dataclass(frozen=True)
class EndToStation:
    kind: Literal["end_to_station"] = "end_to_station"
    end: EndRef
    target: StationRef
    dofs: Optional[Mask6] = None  # None = all relevant DOFs

Connection = Union[EndToStation]  # later add MPC, rigid links, etc.

@dataclass
class Frame:
    nodes: Dict[NodeId, Node] = field(default_factory=dict)
    members: Dict[MemberId, Member] = field(default_factory=dict)
    connections: List[Connection] = field(default_factory=list)
    # Registries
    materials: Dict[str, Any] = field(default_factory=dict) # Replace Any with Material type
    sections: Dict[str, Any] = field(default_factory=dict)  # Replace Any with SectionRegistry type

    # SectionRegistry should implement __getitem__ to auto-load standard shapes (e.g. AISC)
    # if key missing but looks like "W14x90"

    merge_tol: float = 1e-6

    def connect_end_to_member(self, end: EndRef, target: StationRef, dofs: Optional[Mask6] = None) -> None:
        self.connections.append(EndToStation(end=end, target=target, dofs=dofs))
        # build()/analyze() should ensure any referenced StationRef exists in Member.stations.

    def validate(self) -> None:
        pass

    def analyze(self, load_case, settings=None):
        # Local import to enforce one-way dependency (Level 2 -> Level 4)
        from beamy.analysis.engine import run_analysis
        return run_analysis(self, load_case, settings)

    def plot(self, result=None, **kwargs):
        pass
```

(Implementation detail: if you keep `Member` frozen, you’ll update `stations` by replacing the Member in `Frame.members` with a new instance.)

---

# Why this is “better practice”

* Connections are explicit and validateable.
* Meshing points are cleanly separated from loads.
* Second-order is handled in solver state, not polluting your core model.
* Member-only analysis is just a special case of Frame analysis.

If you answer one thing: **2D or 3D first?**
I’ll give you the exact recommended DOF masks, connection DOF options (e.g., pins in 2D), and the minimal `LoadCase` classes that match your solver plan.
