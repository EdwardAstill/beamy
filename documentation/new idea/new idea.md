

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

> I’d rename `node_i/node_j` to `end_i/end_j` to make “member ends” a first-class concept.

### B) Properties (structural)

* `material: MaterialRef`
* `section: SectionRef`
* `orientation: OrientationDef`
  One of:

  * `local_z: Vec3`
  * `roll_angle: float`
  * `ref_node: NodeId`

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

* `geom_nonlinear: bool = True`
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

* `nodes: dict[NodeId, Node]`
* `members: dict[MemberId, Member]`

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

During `Frame.connect_end_to_member(...)`, you:

* add a `Connection(EndToStation(...))`
* also add the station `s` into the target member’s `stations` set (so analysis knows it must split there)

### 2) Remove `Connection.data: dict`

Use typed fields. You’ll thank yourself later when validating and debugging.

### 3) Add second-order readiness to Member

Minimal flags: `geom_nonlinear` + `formulation`.

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

DofMask6 = Tuple[bool, bool, bool, bool, bool, bool]  # True=released (member end) OR restrained (support) depending on context

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

@dataclass(frozen=True)
class Member:
    id: MemberId
    end_i: NodeId
    end_j: NodeId
    kind: Literal["beam", "truss", "cable"]
    material_id: str
    section_id: str
    orientation: Orientation = Orientation()
    releases_i: DofMask6 = (False, False, False, False, False, False)
    releases_j: DofMask6 = (False, False, False, False, False, False)

    # meshing stations (0..1), no loads here
    stations: Tuple[float, ...] = ()

    # second-order controls
    geom_nonlinear: bool = True
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
    end: EndRef = None           # type: ignore
    target: StationRef = None    # type: ignore
    dofs: Optional[DofMask6] = None  # None = all relevant DOFs

Connection = Union[EndToStation]  # later add MPC, rigid links, etc.

@dataclass
class Frame:
    nodes: Dict[NodeId, Node] = field(default_factory=dict)
    members: Dict[MemberId, Member] = field(default_factory=dict)
    connections: List[Connection] = field(default_factory=list)

    merge_tol: float = 1e-6

    def connect_end_to_member(self, end: EndRef, target: StationRef, dofs: Optional[DofMask6] = None) -> None:
        self.connections.append(EndToStation(end=end, target=target, dofs=dofs))
        # also ensure a split station exists on target member (copy-on-write if using frozen dataclass)
        # (implementation detail)

    def validate(self) -> None:
        pass

    def analyze(self, load_case, settings=None):
        pass

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
