# Beamy API Reference

## Overview

Beamy is a structural analysis library for 1D beam and 3D frame analysis. The main public API is organized into several categories:
- **Core Objects**: Materials, supports, loads, and results
- **1D Beam Analysis**: Simple beam models with analysis
- **3D Frame Analysis**: Complex frame structures with nodes and members
- **Visualization**: Plotting utilities for results and geometry

---

## Core Classes

### Material

Represents material properties for beam/member analysis.

```python
@dataclass
class Material:
    name: str                      # Material name (e.g., "Steel")
    E: float                       # Young's modulus
    G: float                       # Shear modulus
    Fy: Optional[float] = None     # Yield stress (for AISC checks)
    transparency: bool = False     # Plot transparency flag
```

**Example:**
```python
from beamy import Material
steel = Material(name="A36 Steel", E=29000, G=11200, Fy=36000)
```

---

### Support

Represents boundary conditions (fixed DOFs) at a location on a beam.

**Key Functions:**
- `validate_support_type(support: str) -> str`: Validates a 6-character support string
- `validate_support_pairs(supports: List[Support]) -> None`: Ensures no duplicate support locations

**Support String Format:**
- 6 characters: `"XYZRXYRZ"` where each is `"0"` (free) or `"1"` (fixed)
- Positions: X, Y, Z (translations), RX, RY, RZ (rotations)

**Example:**
```python
from beamy import Support
pin_support = Support(x=0.0, type="111000")     # Pin support (fixed X, Y, Z)
roller_support = Support(x=10.0, type="011000") # Roller (fixed Y, Z)
```

---

### Load Classes

#### PointForce

Point force applied at a coordinate on a 1D beam.

```python
@dataclass
class PointForce:
    point: np.ndarray      # [x, y, z] coordinates
    force: np.ndarray      # [Fx, Fy, Fz] force vector
```

#### Moment

Concentrated moment applied at position x on a beam.

```python
@dataclass
class Moment:
    x: float               # Position along beam
    moment: np.ndarray     # [Mx, My, Mz] moment vector
```

#### DistributedForce

Distributed load (linear variation) on a 1D beam.

```python
@dataclass
class DistributedForce:
    start_position: np.ndarray     # [x, y, z] start point
    end_position: np.ndarray       # [x, y, z] end point
    start_force: np.ndarray        # [Fx, Fy, Fz] at start
    end_force: np.ndarray          # [Fx, Fy, Fz] at end
```

#### NodalForce

Force applied directly at a node in 3D frames.

```python
@dataclass
class NodalForce:
    node_id: str                           # Target node ID
    force: np.ndarray                      # [FX, FY, FZ] force vector
    coords: str = "global"                 # "global" or "local"
    reference_member_id: Optional[str] = None  # Required if coords="local"
```

#### NodalMoment

Moment applied at a node in 3D frames.

```python
@dataclass
class NodalMoment:
    node_id: str                           # Target node ID
    moment: np.ndarray                     # [MX, MY, MZ] moment vector
    coords: str = "global"                 # "global" or "local"
    reference_member_id: Optional[str] = None  # Required if coords="local"
```

#### MemberPointForce

Point force applied along a frame member.

```python
@dataclass
class MemberPointForce:
    member_id: str                   # Target member ID
    position: float                  # Position along member
    force: np.ndarray                # [FX, FY, FZ] force in member coords
    coords: str = "local"            # "local" or "global"
    position_type: str = "absolute"  # "absolute" (units) or "relative" (0-1)
```

#### MemberDistributedForce

Distributed load over a frame member segment.

```python
@dataclass
class MemberDistributedForce:
    member_id: str           # Target member ID
    start_position: float    # Start position along member
    end_position: float      # End position (-1.0 = full length)
    start_force: np.ndarray  # [FX, FY, FZ] at start
    end_force: np.ndarray    # [FX, FY, FZ] at end
    coords: str = "local"    # "local" or "global"
```

#### NodalSpring

Elastic spring stiffness attached at a node (for elastic soil/foundation).

```python
@dataclass
class NodalSpring:
    node_id: str                           # Target node ID
    K: np.ndarray                          # 6x6 stiffness matrix
    coords: str = "global"                 # "global" or "local"
    reference_member_id: Optional[str] = None  # Required if coords="local"
```

#### LoadCase

Container for loads applied to a 1D beam.

```python
@dataclass
class LoadCase:
    name: str                                  # Load case name
    point_forces: List[PointForce] = []        # Point forces
    moments: List[Moment] = []                 # Concentrated moments
    dist_forces: List[DistributedForce] = []   # Distributed loads
    dist_force_resolution: int = 11            # Mesh points for distrib. loads
    
    def add_point_force(self, pf: PointForce) -> None
    def add_moment(self, m: Moment) -> None
    def add_distributed_force(self, df: DistributedForce) -> None
```

#### FrameLoadCase

Container for loads applied to a 3D frame.

```python
@dataclass
class FrameLoadCase:
    name: str
    nodal_forces: List[NodalForce] = []
    nodal_moments: List[NodalMoment] = []
    member_point_forces: List[MemberPointForce] = []
    member_point_moments: List[MemberPointMoment] = []
    member_dist_forces: List[MemberDistributedForce] = []
    nodal_springs: List[NodalSpring] = []
```

---

### Result / AnalysisResult

Base container for analysis results.

```python
@dataclass
class Result:
    name: str
    values: np.ndarray

@dataclass
class AnalysisResult:
    name: str
    results: Dict[str, Result]
```

---

## 1D Beam Analysis

### Beam1D

Represents a single prismatic beam.

```python
@dataclass
class Beam1D:
    L: float                  # Length
    material: Material        # Material properties
    section: Section          # Cross-section (from sectiony library)
    supports: List[Support]   # Boundary conditions
```

**Example:**
```python
from beamy import Beam1D, Material, Support
from sectiony import Section

steel = Material(name="Steel", E=29000, G=11200, Fy=36000)
section = Section.i_section(height=12, width=6.5, tw=0.3, tf=0.4)
beam = Beam1D(
    L=20.0,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),
        Support(x=20.0, type="011000")
    ]
)
```

---

### LoadedMember

Analysis result for a loaded 1D beam.

```python
@dataclass
class LoadedMember:
    beam: Beam1D
    load_case: LoadCase
    
    # Properties with cached computation:
    @property
    def reactions: Dict[str, float]
    @property
    def support_reactions: List[Tuple[float, Dict[str, float]]]
    @property
    def internal_forces(self) -> Dict[str, List[Tuple[float, float]]]
    @property
    def deflections(self) -> Dict[str, Callable[[float], float]]
    @property
    def stresses(self) -> Dict[str, Callable[[float], float]]
    @property
    def utilization(self) -> float
    @property
    def max_bending_stress: float
    @property
    def max_shear_stress: float
```

**Example:**
```python
from beamy import LoadedMember, PointForce

load_case = LoadCase(name="Dead Load")
load_case.add_point_force(PointForce(point=[10.0, 0, 0], force=[0, -10, 0]))

analysis = LoadedMember(beam=beam, load_case=load_case)
print(f"Max deflection: {analysis.max_deflection}")
print(f"Max stress: {analysis.max_bending_stress}")
```

---

## 3D Frame Analysis

### Node

A point in 3D space where members connect.

```python
@dataclass
class Node:
    id: str                    # Unique node identifier
    position: np.ndarray       # [x, y, z] coordinates (3D)
    support: Optional[str] = None  # Support string (6 chars)
    connected_members: List[str] = []  # Connected member IDs (auto-populated)
    
    # Read-only properties:
    @property
    def __repr__(self) -> str  # String representation with position
```

**Example:**
```python
from beamy import Node
import numpy as np

node_a = Node(id="A", position=np.array([0, 0, 0]), support="111111")
node_b = Node(id="B", position=np.array([10, 0, 0]))
```

---

### Member

A beam element connecting two nodes.

```python
@dataclass
class Member:
    id: str                         # Unique member identifier
    start_node_id: str              # Start node ID
    end_node_id: str                # End node ID
    section: Section                # Cross-section
    material: Material              # Material
    orientation: np.ndarray         # Local Y-axis direction [3D vector]
    element_type: str = "beam"      # "beam", "truss", or "cable"
    releases: Optional[str] = None  # End releases (12-char string)
    constraints: Optional[str] = None  # DOF constraints (12-char string)
    
    # Read-only properties:
    @property
    def length(self) -> float
    @property
    def direction(self) -> np.ndarray  # Normalized X-axis [3D]
    @property
    def local_y(self) -> np.ndarray    # Local Y-axis [3D]
    @property
    def local_z(self) -> np.ndarray    # Local Z-axis [3D]
    @property
    def transformation_matrix(self) -> np.ndarray  # 3x3 rotation matrix
```

**Example:**
```python
from beamy import Member
import numpy as np

member = Member(
    id="M1",
    start_node_id="A",
    end_node_id="B",
    section=section,
    material=steel,
    orientation=np.array([0, 1, 0])  # Local Y points in global Y
)
```

---

### Frame

A collection of nodes and members forming a 3D structure.

```python
@dataclass
class Frame:
    members: List[Member]
    nodes: Dict[str, Node]  # Auto-populated from members
    
    @classmethod
    def from_nodes_and_members(cls, nodes: List[Node], members: List[Member]) -> Frame
        """Factory method to create a Frame."""
    
    # Properties:
    @property
    def node_positions(self) -> Dict[str, np.ndarray]
    @property
    def supported_nodes(self) -> List[Node]
    @property
    def member_lengths(self) -> Dict[str, float]
    
    # Methods:
    def get_member(self, mid: str) -> Member
    def get_node(self, nid: str) -> Node
    def members_at_node(self, nid: str) -> List[Member]
```

**Example:**
```python
from beamy import Frame

nodes = [node_a, node_b, node_c]
members = [member1, member2]
frame = Frame.from_nodes_and_members(nodes, members)

member = frame.get_member("M1")
print(f"Member length: {member.length}")
```

---

### FrameAnalysisSettings

Configuration for frame analysis method and convergence controls.

```python
@dataclass
class FrameAnalysisSettings:
    analysis_method: str = "GENERIC_LINEAR"
        # Options: "GENERIC_LINEAR", "SECOND_ORDER_ELASTIC", 
        #          "AISC360_DAM", "EC3_GLOBAL_IMPERFECTIONS", "AS4100_SECOND_ORDER"
    
    # Second-order effects:
    imperfection_model: str = "none"       # "none", "notional_loads", "initial_sway"
    notional_factor: float = 0.002
    notional_p_source: str = "input_based" # "input_based" or "reactions_based"
    notional_axes: Tuple = ("x", "y")
    notional_signs: Tuple[float, float] = (1.0, 1.0)
    
    # Stiffness reduction rules:
    stiffness_rules: str = "none"  # "none", "aisc_dam", "ec3", "as4100"
    bending_stiffness_factor: Optional[float] = None
    torsion_stiffness_factor: float = 1.0
    axial_stiffness_factor: float = 1.0
    
    # Convergence:
    max_iter: int = 25
    tol_u_rel: float = 1e-6
    tol_u_abs: float = 1e-12
    n_steps: int = 1
    relaxation_omega: float = 1.0
    
    # Element behavior:
    cable_tension_only: bool = True
```

---

### LoadedFrame

Analysis result for a loaded frame structure.

```python
@dataclass
class LoadedFrame:
    frame: Frame
    load_case: FrameLoadCase
    settings: Optional[FrameAnalysisSettings] = None
    
    # Execute analysis:
    def analyze(self) -> FrameAnalysisResult
    
    # Convenience properties (call analyze() internally):
    @property
    def nodal_displacements(self) -> Dict[str, np.ndarray]
    @property
    def reactions(self) -> Dict[str, np.ndarray]
    @property
    def member_end_forces(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]
```

**Example:**
```python
from beamy import LoadedFrame, FrameLoadCase, NodalForce

load_case = FrameLoadCase(name="Live Load")
load_case.nodal_forces.append(NodalForce(node_id="B", force=np.array([0, -5, 0])))

analysis = LoadedFrame(frame=frame, load_case=load_case)
result = analysis.analyze()

print(f"Node B displacement: {result.nodal_displacements['B']}")
print(f"Reactions: {result.reactions}")
```

---

### FrameAnalysisResult

Primary output container from frame analysis.

```python
@dataclass
class FrameAnalysisResult:
    nodal_displacements: Dict[str, np.ndarray]      # [Ux, Uy, Uz, Rx, Ry, Rz]
    reactions: Dict[str, np.ndarray]                # [Rx, Ry, Rz, Mx, My, Mz]
    reactions_physical: Dict[str, np.ndarray]       # Physical support reactions
    reactions_stabilization: Dict[str, np.ndarray]  # Numerical stabilization only
    equilibrium_residual_norm: float
    member_end_forces: Dict[str, Tuple[np.ndarray, np.ndarray]]  # Start and end forces
    demand_provider: Optional[MemberDemandProvider]
    
    # Metadata:
    settings: FrameAnalysisSettings
    converged: bool
    iterations: int
    design_grade_ok: bool
```

---

### MemberResults

Results for a single frame member after analysis.

```python
@dataclass
class MemberResults:
    member_id: str
    start_forces: np.ndarray      # [Fx, Fy, Fz, Mx, My, Mz] at start node
    end_forces: np.ndarray        # [Fx, Fy, Fz, Mx, My, Mz] at end node
    
    # Convenience properties:
    @property
    def start_axial_force(self) -> float
    @property
    def start_shear_y(self) -> float
    @property
    def start_shear_z(self) -> float
    @property
    def start_torsion(self) -> float
    @property
    def start_moment_y(self) -> float
    @property
    def start_moment_z(self) -> float
    # ... similar for end_* properties
```

---

## Visualization Functions

### Beam Plotting

#### plot_beam_diagram

Plots a 1D beam with supports and loads.

```python
def plot_beam_diagram(
    beam: Beam1D,
    load_case: Optional[LoadCase] = None,
    show_section: bool = True,
    figsize: Tuple[float, float] = (14, 6)
) -> None
```

#### plot_analysis_results

Plots beam internal forces and deflections.

```python
def plot_analysis_results(
    loaded_beam: LoadedMember,
    figsize: Tuple[float, float] = (14, 10)
) -> None
```

### Frame Plotting

#### plot_frame

Plots a 3D frame structure with nodes, members, and supports.

```python
def plot_frame(
    frame: Frame,
    show_ids: bool = True,
    figsize: Tuple[float, float] = (10, 8)
) -> None
```

#### plot_deflection

Plots deformed frame geometry after analysis.

```python
def plot_deflection(
    loaded_frame: LoadedFrame,
    scale: float = 1.0,
    show_original: bool = True
) -> None
```

#### plot_von_mises

Plots Von Mises stress distribution on frame members.

```python
def plot_von_mises(
    loaded_frame: LoadedFrame,
    scale: float = 1.0
) -> None
```

#### plot_member_diagrams

Plots internal force/moment diagrams for frame members.

```python
def plot_member_diagrams(
    loaded_frame: LoadedFrame,
    member_ids: Optional[List[str]] = None
) -> None
```

### Section & Support Plotting

#### plot_section

Visualizes a cross-section geometry.

```python
def plot_section(section: Section) -> None
```

#### plot_supports

Shows support conditions at nodes/members.

```python
def plot_supports(frame: Frame) -> None
```

#### plot_loads

Visualizes applied loads on a structure.

```python
def plot_loads(
    structure,  # Beam1D or Frame
    load_case   # LoadCase or FrameLoadCase
) -> None
```

### StressPlotter

Advanced stress visualization utility.

```python
class StressPlotter:
    def __init__(self, loaded_beam: LoadedMember)
    
    def plot_bending_stress(self, x: float) -> None
        """Plot bending stress distribution at cross-section at position x."""
    
    def plot_shear_stress(self, x: float) -> None
        """Plot shear stress distribution at position x."""
    
    def plot_von_mises(self) -> None
        """Plot Von Mises stress envelope along beam."""
```

---

## Utility Functions

### Frame Builder

#### FrameBuilder

Builder pattern for constructing frames programmatically.

```python
class FrameBuilder:
    def add_node(self, id: str, x: float, y: float, z: float, support: Optional[str] = None) -> "FrameBuilder"
    def add_member(self, id: str, start_id: str, end_id: str, section: Section, material: Material, orientation: np.ndarray) -> "FrameBuilder"
    def build(self) -> Frame
```

**Example:**
```python
from beamy import FrameBuilder

builder = FrameBuilder()
builder.add_node("A", 0, 0, 0, support="111111")
builder.add_node("B", 10, 0, 0)
builder.add_node("C", 10, 10, 0, support="111111")
builder.add_member("M1", "A", "B", section, steel, np.array([0, 1, 0]))
builder.add_member("M2", "B", "C", section, steel, np.array([-1, 0, 0]))
frame = builder.build()
```

#### round_coord

Utility to round coordinates for node merging.

```python
def round_coord(value: float, decimals: int = 6) -> float
```

---

## Checking (AISC 9)

### aisc_9 Module

AISC 360-16 strength checks (Chapter 9).

```python
def check_member(
    member: Member,
    analysis_result: FrameAnalysisResult,
    phi_b: float = 0.9
) -> Dict[str, float]
    """Return utilization ratios for each check."""
```

---

## External Dependencies

**sectiony**: Provides `Section` and `Geometry` classes for cross-section properties.

```python
from sectiony import Section, Geometry

# Create sections:
I_section = Section.i_section(height=12, width=6.5, tw=0.3, tf=0.4)
rect_section = Section.rectangular(width=6, height=12)
circle_section = Section.circular(diameter=4)
```

---

## Common Workflows

### Analyze a Simple Cantilever Beam

```python
from beamy import *
from sectiony import Section
import numpy as np

# Setup
steel = Material(name="Steel", E=29000, G=11200, Fy=36000)
section = Section.i_section(height=12, width=6.5, tw=0.3, tf=0.4)

beam = Beam1D(
    L=20,
    material=steel,
    section=section,
    supports=[Support(x=0, type="111111")]
)

# Load
loads = LoadCase(name="Gravity")
loads.add_point_force(PointForce(point=[20, 0, 0], force=[0, -10, 0]))

# Analyze
analysis = LoadedMember(beam=beam, load_case=loads)
print(f"Max deflection: {max(analysis.deflections['z']())}")
print(f"Max stress: {analysis.max_bending_stress}")

# Plot
plot_analysis_results(analysis)
```

### Analyze a 2D Portal Frame

```python
# Build frame
builder = FrameBuilder()
builder.add_node("A", 0, 0, 0, support="111111")
builder.add_node("B", 0, 20, 0)
builder.add_node("C", 30, 20, 0)
builder.add_node("D", 30, 0, 0, support="111111")

builder.add_member("M1", "A", "B", section, steel, np.array([0, 1, 0]))
builder.add_member("M2", "B", "C", section, steel, np.array([1, 0, 0]))
builder.add_member("M3", "C", "D", section, steel, np.array([0, -1, 0]))

frame = builder.build()

# Load
load_case = FrameLoadCase(name="Wind")
load_case.nodal_forces.append(NodalForce(node_id="B", force=np.array([10, 0, 0])))
load_case.nodal_forces.append(NodalForce(node_id="C", force=np.array([10, 0, 0])))

# Analyze
settings = FrameAnalysisSettings(analysis_method="GENERIC_LINEAR")
analysis = LoadedFrame(frame=frame, load_case=load_case, settings=settings)
result = analysis.analyze()

print(f"Reactions: {result.reactions}")

# Plot
plot_frame(frame)
plot_deflection(analysis)
```
