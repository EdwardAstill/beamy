## Beamy API Reference (Proposed vNext)

**Status**: This document defines the **proposed** public API for a future Beamy release. It is a design draft and **may not match the current implementation** in `src/beamy` yet.

### Design goals

- **One load container**: `LoadCase` is used everywhere.
- **No separate 1D beam type**: remove `Beam1D`.
- **Single-member analysis is first-class**: `LoadedMember` analyzes a member in isolation (internally equivalent to a 2-node / 1-member structure).
- **Frame analysis is explicit**: users call `frame.analyze(load_case, settings=...)` and results are cached on the `Frame`.
- **Applied loads vs demand**: a member can have non-trivial internal forces under a `LoadCase` **even if the load case contains no direct loads on that member**, due to global compatibility and load-path effects.

---

### Quick start

#### Analyze one member “in isolation”

```python
from __future__ import annotations

import numpy as np
from sectiony import Section

from beamy import (
    Material,
    LoadCase,
    MemberPointForce,
    LoadedMember,
)

steel = Material(name="A36", E=29000.0, G=11200.0, Fy=36000.0)
section = Section.i_section(height=12.0, width=6.5, tw=0.3, tf=0.4)

load_case = LoadCase(name="LC1")
load_case.member_point_forces.append(
    MemberPointForce(
        member_id="M1",
        position=10.0,
        force=np.array([0.0, -10.0, 0.0]),
        coords="global",
        position_type="absolute",
    )
)

beam = LoadedMember(
    id="M1",
    start=np.array([0.0, 0.0, 0.0]),
    end=np.array([20.0, 0.0, 0.0]),
    section=section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0]),
    support_start="111000",
    support_end="011000",
    load_case=load_case,
)

result = beam.analyze()
profile = beam.member_demand().actions(points=201)
```

#### Analyze a frame (multi-member)

```python
from __future__ import annotations

import numpy as np
from sectiony import Section

from beamy import (
    Material,
    Node,
    Member,
    Frame,
    LoadCase,
    NodalForce,
)

steel = Material(name="Steel", E=29000.0, G=11200.0, Fy=36000.0)
section = Section.i_section(height=12.0, width=6.5, tw=0.3, tf=0.4)

members = [
    Member(id="M1", start=np.array([0.0, 0.0, 0.0]), end=np.array([10.0, 0.0, 0.0]), section=section, material=steel, orientation=np.array([0.0, 0.0, 1.0])),
    Member(id="M2", start=np.array([10.0, 0.0, 0.0]), end=np.array([10.0, 10.0, 0.0]), section=section, material=steel, orientation=np.array([0.0, 0.0, 1.0])),
]

frame = Frame.from_members(members)

# Apply supports to auto-generated nodes
for node in frame.nodes.values():
    if np.allclose(node.position, [0.0, 0.0, 0.0]):
        node.support = "111111"
    elif np.allclose(node.position, [10.0, 10.0, 0.0]):
        node.support = "111111"

load_case = LoadCase(name="Wind")
load_case.nodal_forces.append(NodalForce(node_id="B", force=np.array([10.0, 0.0, 0.0]), coords="global"))

result = frame.analyze(load_case)

# Members can have demand even if they were not directly loaded in the load_case.
m2 = frame.member_demand("M2")
envelopes = m2.envelopes(points=201)
```

---

### Core modeling

#### `Material`

```python
from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class Material:
    name: str
    E: float
    G: float
    Fy: Optional[float] = None
```

#### Support strings (boundary conditions)

Support codes are 6-character `"0"`/`"1"` strings representing fixed/free DOFs:

`[UX, UY, UZ, RX, RY, RZ]`

Examples:
- `"111111"` fixed
- `"111000"` pinned (translations fixed, rotations free)
- `"011000"` roller (UY/UZ fixed)

#### `Node` (frame)

```python
from dataclasses import dataclass, field
from typing import Optional, List
import numpy as np

@dataclass
class Node:
    id: str
    position: np.ndarray  # (3,)
    support: Optional[str] = None
    connected_members: List[str] = field(default_factory=list)
```

#### `Member`

```python
from dataclasses import dataclass
from typing import Optional, Literal
import numpy as np
from sectiony import Section

@dataclass
class Member:
    id: str
    start: np.ndarray           # (3,) start position [x, y, z]
    end: np.ndarray             # (3,) end position [x, y, z]
    section: Section
    material: Material
    orientation: np.ndarray     # local Y direction, (3,)
    element_type: Literal["beam", "truss", "cable"] = "beam"
    releases: Optional[str] = None      # 12-digit release string
    constraints: Optional[str] = None   # 12-digit constraint string
```

**Note**: Connectivity is automatic. The `Frame` generates nodes from member endpoints and merges coincident positions.

#### `Frame`

```python
from dataclasses import dataclass
from typing import Dict, List

@dataclass
class Frame:
    members: List[Member]
    nodes: Dict[str, Node]  # auto-generated

    @classmethod
    def from_members(cls, members: List[Member], node_tolerance: float = 1e-6) -> "Frame":
        """Build frame from members, auto-generating nodes at endpoints."""
        ...

    @classmethod
    def from_nodes_and_members(cls, nodes: List[Node], members: List[Member]) -> "Frame":
        """Build frame from explicit nodes and members (legacy/advanced)."""
        ...
```

**Recommended**: Use `Frame.from_members(members)` - it auto-generates nodes and connectivity.

---

### Loads

#### `LoadCase` (used for both standalone members and frames)

`LoadCase` stores **applied loads**. After analysis, you query **member demand** (actions/stresses/deflections) for any member.

```python
from dataclasses import dataclass, field
from typing import List

@dataclass
class LoadCase:
    name: str

    nodal_forces: List["NodalForce"] = field(default_factory=list)
    nodal_moments: List["NodalMoment"] = field(default_factory=list)
    nodal_springs: List["NodalSpring"] = field(default_factory=list)

    member_point_forces: List["MemberPointForce"] = field(default_factory=list)
    member_point_moments: List["MemberPointMoment"] = field(default_factory=list)
    member_distributed_forces: List["MemberDistributedForce"] = field(default_factory=list)

    # “Supports as loads”: useful for intermediate restraints without manual member splitting
    member_point_supports: List["MemberPointSupport"] = field(default_factory=list)
    member_supports: List["MemberSupport"] = field(default_factory=list)
```

#### Nodal loads

```python
from dataclasses import dataclass
from typing import Optional, Literal
import numpy as np

@dataclass(frozen=True)
class NodalForce:
    node_id: str
    force: np.ndarray  # (3,) [FX, FY, FZ]
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None  # required if coords="local"

@dataclass(frozen=True)
class NodalMoment:
    node_id: str
    moment: np.ndarray  # (3,) [MX, MY, MZ]
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None  # required if coords="local"
```

#### Member loads

```python
from dataclasses import dataclass
from typing import Literal, Optional
import numpy as np

@dataclass(frozen=True)
class MemberPointForce:
    member_id: str
    position: float
    force: np.ndarray  # (3,) [FX, FY, FZ]
    coords: Literal["global", "local"] = "local"
    position_type: Literal["absolute", "relative"] = "absolute"

@dataclass(frozen=True)
class MemberPointMoment:
    member_id: str
    position: float
    moment: np.ndarray  # (3,) [MX, MY, MZ]
    coords: Literal["global", "local"] = "local"
    position_type: Literal["absolute", "relative"] = "absolute"

@dataclass(frozen=True)
class MemberDistributedForce:
    member_id: str
    start_position: float
    end_position: Optional[float]  # None means “to member end”
    start_force: np.ndarray        # (3,) per length
    end_force: np.ndarray          # (3,) per length
    coords: Literal["global", "local"] = "local"
```

#### Springs (optional)

```python
from dataclasses import dataclass
from typing import Optional, Literal
import numpy as np

@dataclass(frozen=True)
class NodalSpring:
    node_id: str
    K: np.ndarray  # (6,6) in [UX, UY, UZ, RX, RY, RZ]
    coords: Literal["global", "local"] = "global"
    reference_member_id: Optional[str] = None  # required if coords="local"
```

#### Member supports via `LoadCase` (optional)

These are convenient when you want to add restraints along a member (or on all nodes along a member) without manually splitting members.

```python
from dataclasses import dataclass
from typing import Literal

@dataclass(frozen=True)
class MemberPointSupport:
    member_id: str
    position: float
    support: str  # 6-digit support code
    position_type: Literal["absolute", "relative"] = "absolute"

@dataclass(frozen=True)
class MemberSupport:
    member_id: str
    support: str  # 6-digit support code
```

---

### Analysis

#### `FrameAnalysisSettings` (shared by member + frame analysis)

```python
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

AnalysisMethod = Literal[
    "GENERIC_LINEAR",
    "SECOND_ORDER_ELASTIC",
    "AISC360_DAM",
    "EC3_GLOBAL_IMPERFECTIONS",
    "AS4100_SECOND_ORDER",
]

@dataclass
class FrameAnalysisSettings:
    analysis_method: AnalysisMethod = "GENERIC_LINEAR"
    imperfection_model: Literal["none", "notional_loads", "initial_sway"] = "none"
    notional_factor: float = 0.002
    notional_axes: Tuple[Literal["x", "y"], ...] = ("x", "y")

    stiffness_rules: Literal["none", "aisc_dam", "ec3", "as4100"] = "none"
    bending_stiffness_factor: Optional[float] = None
    torsion_stiffness_factor: float = 1.0
    axial_stiffness_factor: float = 1.0

    max_iter: int = 25
    tol_u_rel: float = 1e-6
    tol_u_abs: float = 1e-12
    n_steps: int = 1
    relaxation_omega: float = 1.0

    cable_tension_only: bool = True
```

#### `LoadedMember` (single-member analysis)

`LoadedMember` is the “beam” object in vNext: it is a **member + end conditions + load_case** that can be analyzed in isolation.

```python
import numpy as np
from sectiony import Section
from typing import Optional

class LoadedMember:
    def __init__(
        self,
        *,
        id: str,
        start: np.ndarray,
        end: np.ndarray,
        section: Section,
        material: Material,
        orientation: np.ndarray,
        support_start: str,
        support_end: str,
        load_case: LoadCase,
        settings: Optional[FrameAnalysisSettings] = None,
    ) -> None: ...

    def analyze(self) -> "FrameAnalysisResult": ...
    def member_demand(self) -> "MemberDemand": ...
```

#### `Frame.analyze(...)` (multi-member analysis)

```python
from typing import Optional

class Frame:
    def analyze(self, load_case: LoadCase, settings: Optional[FrameAnalysisSettings] = None) -> "FrameAnalysisResult": ...
    def member_demand(self, member_id: str) -> "MemberDemand": ...
```

---

### Results & demand

#### `Result` and `MemberActionProfile`

```python
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Result:
    x: np.ndarray
    values: np.ndarray

    def at(self, x_loc: float) -> float: ...
    @property
    def max(self) -> float: ...
    @property
    def min(self) -> float: ...
    @property
    def abs_max(self) -> float: ...

@dataclass(frozen=True)
class MemberActionProfile:
    member_id: str
    length: float

    axial: Result
    shear_y: Result
    shear_z: Result
    torsion: Result
    bending_y: Result
    bending_z: Result
```

#### `MemberDemand`

`MemberDemand` represents the **solved** response of a member under the structure’s `load_case`.

```python
from typing import Dict, Tuple

class MemberDemand:
    def actions(self, points: int = 201) -> MemberActionProfile: ...
    def envelopes(self, points: int = 201) -> Dict[str, Tuple[float, float]]: ...
```

#### `FrameAnalysisResult`

`FrameAnalysisResult` is the primary output container from analysis.

```python
from dataclasses import dataclass
from typing import Dict, Optional, Tuple
import numpy as np

@dataclass(frozen=True)
class FrameAnalysisResult:
    nodal_displacements: Dict[str, np.ndarray]  # [UX, UY, UZ, RX, RY, RZ]
    reactions: Dict[str, np.ndarray]            # [FX, FY, FZ, MX, MY, MZ]
    member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]]  # local start/end forces
    demand_provider: Optional[object]
```

---

### Plotting (proposed signatures)

All plotting should be consistent with the analysis objects and should support saving to **`.svg`**.

```python
from typing import Optional

def plot_frame(frame: Frame, load_case: Optional[LoadCase] = None, *, save_path: Optional[str] = None) -> None: ...
def plot_deflection(frame: Frame, *, scale_factor: float = 1.0, save_path: Optional[str] = None) -> None: ...
def plot_von_mises(frame: Frame, *, save_path: Optional[str] = None) -> None: ...
def plot_member_diagrams(frame: Frame, member_id: str, *, save_path: Optional[str] = None) -> None: ...
```

Example:

```python
frame.analyze(load_case)
plot_deflection(frame, scale_factor=25.0, save_path="deflection.svg")
```

---

### Checks (AISC 360 Chapter 9)

Checks should run on **solved member demand**, not on the raw loads.

Proposed interface:

```python
def aisc_9_check(profile: MemberActionProfile, *, length_unit: str, force_unit: str): ...
```

Example (frame member):

```python
frame.analyze(load_case)

demand = frame.member_demand("M7")
profile = demand.actions(points=801)

check = aisc_9_check(profile, length_unit="m", force_unit="N")
utilization = check.utilisation
```

---

### Migration from current API (high-level)

- **`Beam1D` removed**:
  - Old: `Beam1D(...)` + `LoadedMember(beam=..., loads=...)`
  - New: `LoadedMember(..., load_case=...)`
- **`FrameLoadCase` → `LoadCase`**:
  - Old: `FrameLoadCase(name="...")`
  - New: `LoadCase(name="...")`
- **`LoadedFrame` / `FrameAnalysis` removed from the public API**:
  - New: call `frame.analyze(load_case=..., settings=...)` and query demand via `frame.member_demand(...)`
- **Use `load_case` everywhere**:
  - Replace any `loads=` usage with `load_case=` (argument name + attribute name)


