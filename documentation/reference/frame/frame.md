# Frame Module

Reference documentation for 3D frame structures composed of interconnected beam members.

## Table of Contents

- [Overview](#overview)
- [Coordinate Systems](#coordinate-systems)
- [Node](#node)
- [Member](#member)
- [Frame](#frame)
- [FrameLoadCase](#frameloadcase)
- [LoadedFrame](#loadedframe)
  - [Properties](#properties)
  - [Analysis Methods](#analysis-methods)
  - [Visualization](#visualization)
- [Workflow Example](#workflow-example)
- [Design Decisions](#design-decisions)

---

## Overview

A **Frame** represents a 3D structural system composed of multiple beam members connected at nodes. Unlike `Beam1D` which operates along a single axis, a Frame allows members to be oriented arbitrarily in 3D space.

**Key Concepts:**
- **Nodes**: Connection points in 3D space where members meet
- **Members**: Individual beam elements spanning between two nodes
- **Global Coordinate System**: Fixed XYZ reference frame for the entire structure
- **Local Coordinate System**: Each member has its own local xyz axes for section orientation

**Relationship to Existing Code:**
- Reuses `Material` and `Section` from `beamy.setup`
- Internally converts to `LoadedBeam` objects for analysis
- Leverages existing FEM solvers once member forces/moments are resolved

---

## Coordinate Systems

### Global Coordinates (X, Y, Z)

The frame uses a right-handed global coordinate system:
- **X**: Typically horizontal (east/right)
- **Y**: Typically horizontal (north/forward)  
- **Z**: Typically vertical (up)

All node positions are specified in global coordinates.

### Local Member Coordinates (x, y, z)

Each member has a local coordinate system:
- **x**: Along the member axis, from start node to end node
- **y**: Perpendicular to x, defined by the `orientation` vector
- **z**: Perpendicular to both x and y (right-hand rule: z = x × y)

The local y-axis determines how the section is oriented (e.g., which way the web of an I-beam faces).

---

## Node

A point in 3D space where members connect and where supports/loads can be applied.

```python
from beamy.frame import Node

@dataclass
class Node:
    id: str                           # Unique identifier
    position: np.ndarray              # [X, Y, Z] in global coordinates
    support: Optional[str] = None     # 6-digit support string (same format as Support.type)
    connected_members: List[str] = field(default_factory=list)  # Member IDs (auto-populated)
```

### Parameters

- `id` (str): Unique node identifier (e.g., "A", "N1", "base_left")
- `position` (np.ndarray): 3D position vector [X, Y, Z] in global coordinates
- `support` (str, optional): Support constraint string using same format as `Support.type`:
  - Position 1: UX (translation in global X)
  - Position 2: UY (translation in global Y)
  - Position 3: UZ (translation in global Z)
  - Position 4: RX (rotation about global X)
  - Position 5: RY (rotation about global Y)
  - Position 6: RZ (rotation about global Z)
- `connected_members` (List[str]): Automatically populated list of member IDs connected to this node

### Common Support Types

- `"111111"`: Fully fixed (all DOFs constrained)
- `"111000"`: Pinned (translations fixed, rotations free)
- `"110000"`: Roller on Z-plane (UX, UY fixed; UZ and rotations free)
- `"001000"`: Vertical roller (only UZ fixed)
- `None`: No support (internal node or free end)

### Example

```python
import numpy as np
from beamy.frame import Node

# Fixed base support
base = Node(id="A", position=np.array([0.0, 0.0, 0.0]), support="111111")

# Pinned connection at top of column
top = Node(id="B", position=np.array([0.0, 0.0, 5.0]), support=None)

# Roller support
roller = Node(id="C", position=np.array([8.0, 0.0, 5.0]), support="110000")
```

---

## Member

A beam element connecting two nodes in 3D space.

```python
from beamy.frame import Member

@dataclass
class Member:
    id: str                     # Unique identifier
    start_node_id: str          # ID of start node
    end_node_id: str            # ID of end node  
    section: Section            # Cross-section (from sectiony)
    material: Material          # Material properties
    orientation: np.ndarray     # Vector defining local y-axis direction
    releases: Optional[str] = None  # End releases (12-digit string)
```

### Parameters

- `id` (str): Unique member identifier (e.g., "M1", "column_1", "beam_AB")
- `start_node_id` (str): ID of the node at the start of the member
- `end_node_id` (str): ID of the node at the end of the member
- `section` (Section): Cross-section properties (from sectiony library)
- `material` (Material): Material properties (E, G, Fy)
- `orientation` (np.ndarray): A 3D vector in global coordinates that defines the direction of the local y-axis. This vector is projected onto the plane perpendicular to the member axis to determine the actual local y-direction.
- `releases` (str, optional): Member end releases as a 12-digit string (see [End Releases](#end-releases))

### Computed Properties

The following properties are computed automatically:

```python
@property
def length(self) -> float:
    """Member length computed from node positions."""
    
@property
def direction(self) -> np.ndarray:
    """Unit vector from start to end node (local x-axis)."""
    
@property
def local_y(self) -> np.ndarray:
    """Unit vector for local y-axis (from orientation)."""
    
@property
def local_z(self) -> np.ndarray:
    """Unit vector for local z-axis (x × y)."""
    
@property
def transformation_matrix(self) -> np.ndarray:
    """3x3 rotation matrix from local to global coordinates."""
```

### End Releases

The `releases` parameter is a 12-digit string controlling moment/force transfer at member ends:
- Digits 1-6: Start node releases [Nx, Vy, Vz, Tx, My, Mz]
- Digits 7-12: End node releases [Nx, Vy, Vz, Tx, My, Mz]

Where `0` = connected (force/moment transfers) and `1` = released (force/moment does not transfer).

**Common Release Patterns:**
- `None` or `"000000000000"`: Fully connected at both ends (rigid frame behavior)
- `"000001000001"`: Pinned at both ends about z-axis (simple beam)
- `"000001000000"`: Pinned at start only (moment release)

### Orientation Examples

```python
# Vertical column - local y points in global +X direction
column = Member(
    id="col1",
    start_node_id="A", 
    end_node_id="B",
    section=w_section,
    material=steel,
    orientation=np.array([1.0, 0.0, 0.0])  # y-axis → global X
)

# Horizontal beam in XZ plane - local y points up (global +Z)
beam = Member(
    id="beam1",
    start_node_id="B",
    end_node_id="C", 
    section=w_section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0])  # y-axis → global Z
)

# Inclined member - y-axis perpendicular to member, pointing "up"
brace = Member(
    id="brace1",
    start_node_id="A",
    end_node_id="C",
    section=angle_section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0])  # Projected onto perpendicular plane
)
```

---

## Frame

A collection of nodes and members forming a 3D structural system.

```python
from beamy.frame import Frame

@dataclass
class Frame:
    members: List[Member]
    nodes: Dict[str, Node] = field(init=False)  # Auto-built from members
```

### Parameters

- `members` (List[Member]): List of Member objects defining the frame geometry

### Computed Attributes

- `nodes` (Dict[str, Node]): Dictionary mapping node IDs to Node objects. Automatically constructed by:
  1. Extracting unique node IDs from all members
  2. Requiring nodes to be pre-defined OR inferring positions (TBD)
  3. Populating `connected_members` lists

### Alternative Constructor

```python
@classmethod
def from_nodes_and_members(
    cls,
    nodes: List[Node],
    members: List[Member]
) -> Frame:
    """
    Create a Frame from explicit node and member lists.
    
    Validates that all member node references exist in the nodes list.
    """
```

### Validation

On construction, Frame validates:
- All member `start_node_id` and `end_node_id` reference existing nodes
- No duplicate member IDs
- No duplicate node IDs  
- At least one support exists (at least one node has a non-None support)
- Structure has sufficient supports for stability (6 DOFs minimum constrained)
- No zero-length members (start and end nodes at same position)

### Properties and Methods

```python
@property
def node_positions(self) -> Dict[str, np.ndarray]:
    """Dictionary of node ID → position array."""

@property  
def supported_nodes(self) -> List[Node]:
    """List of nodes that have supports defined."""

@property
def member_lengths(self) -> Dict[str, float]:
    """Dictionary of member ID → length."""

def get_member(self, member_id: str) -> Member:
    """Retrieve a member by ID."""

def get_node(self, node_id: str) -> Node:
    """Retrieve a node by ID."""

def members_at_node(self, node_id: str) -> List[Member]:
    """Get all members connected to a node."""
```

### Example

```python
import numpy as np
from sectiony.library import i as i_section
from beamy import Material
from beamy.frame import Node, Member, Frame

# Materials and sections
steel = Material(name="Steel", E=200e9, G=80e9, Fy=250e6)
column_section = i_section(d=0.3, b=0.3, tf=0.02, tw=0.015, r=0.0)
beam_section = i_section(d=0.4, b=0.2, tf=0.015, tw=0.010, r=0.0)

# Define nodes
nodes = [
    Node(id="A", position=np.array([0.0, 0.0, 0.0]), support="111111"),
    Node(id="B", position=np.array([6.0, 0.0, 0.0]), support="111111"),
    Node(id="C", position=np.array([0.0, 0.0, 4.0])),
    Node(id="D", position=np.array([6.0, 0.0, 4.0])),
]

# Define members
members = [
    # Columns
    Member(id="col1", start_node_id="A", end_node_id="C",
           section=column_section, material=steel,
           orientation=np.array([1.0, 0.0, 0.0])),
    Member(id="col2", start_node_id="B", end_node_id="D",
           section=column_section, material=steel,
           orientation=np.array([1.0, 0.0, 0.0])),
    # Beam
    Member(id="beam1", start_node_id="C", end_node_id="D",
           section=beam_section, material=steel,
           orientation=np.array([0.0, 0.0, 1.0])),
]

# Create frame
frame = Frame.from_nodes_and_members(nodes, members)
```

---

## FrameLoadCase

Loads applied to a frame structure.

```python
from beamy.frame import FrameLoadCase

@dataclass
class FrameLoadCase:
    name: str
    nodal_forces: List[NodalForce] = field(default_factory=list)
    nodal_moments: List[NodalMoment] = field(default_factory=list)
    member_point_forces: List[MemberPointForce] = field(default_factory=list)
    member_distributed_forces: List[MemberDistributedForce] = field(default_factory=list)
```

### Load Types

#### NodalForce

Force applied directly at a node (in global coordinates).

```python
@dataclass
class NodalForce:
    node_id: str              # Target node ID
    force: np.ndarray         # [FX, FY, FZ] in global coordinates
```

#### NodalMoment

Moment applied directly at a node (in global coordinates).

```python
@dataclass  
class NodalMoment:
    node_id: str              # Target node ID
    moment: np.ndarray        # [MX, MY, MZ] in global coordinates
```

#### MemberPointForce

Point force applied along a member (in local or global coordinates).

```python
@dataclass
class MemberPointForce:
    member_id: str            # Target member ID
    position: float           # Distance from start node (0 to L) OR fraction (0 to 1)
    force: np.ndarray         # [Fx, Fy, Fz] force vector
    coords: str = "local"     # "local" or "global" coordinate system
    position_type: str = "absolute"  # "absolute" (distance) or "relative" (fraction)
```

#### MemberDistributedForce

Distributed force along a member segment.

```python
@dataclass
class MemberDistributedForce:
    member_id: str
    start_position: float     # Distance from start node
    end_position: float       # Distance from start node  
    start_force: np.ndarray   # [wx, wy, wz] per unit length at start
    end_force: np.ndarray     # [wx, wy, wz] per unit length at end
    coords: str = "local"     # "local" or "global"
```

### Methods

```python
def add_nodal_force(self, node_id: str, force: np.ndarray) -> None:
    """Add a force at a node."""

def add_nodal_moment(self, node_id: str, moment: np.ndarray) -> None:
    """Add a moment at a node."""

def add_member_point_force(
    self,
    member_id: str,
    position: float,
    force: np.ndarray,
    coords: str = "local"
) -> None:
    """Add a point force along a member."""

def add_member_distributed_force(
    self,
    member_id: str,
    start_position: float,
    end_position: float,
    start_force: np.ndarray,
    end_force: np.ndarray,
    coords: str = "local"
) -> None:
    """Add a distributed force along a member segment."""

def add_member_uniform_force(
    self,
    member_id: str,
    force: np.ndarray,
    coords: str = "local"
) -> None:
    """Add a uniform distributed force over the entire member length."""
```

### Example

```python
from beamy.frame import FrameLoadCase
import numpy as np

loads = FrameLoadCase(name="Wind + Dead Load")

# Point load at node D (global coordinates)
loads.add_nodal_force("D", force=np.array([50000.0, 0.0, 0.0]))  # 50kN horizontal

# Uniform dead load on beam (local z = down relative to beam)
loads.add_member_uniform_force("beam1", force=np.array([0.0, 0.0, -5000.0]))  # 5kN/m
```

---

## LoadedFrame

A frame with loads applied, ready for analysis.

```python
from beamy.frame import LoadedFrame

@dataclass
class LoadedFrame:
    frame: Frame
    loads: FrameLoadCase
```

### Initialization

On instantiation, `LoadedFrame`:
1. Validates load references (node IDs and member IDs exist)
2. Transforms member loads from local to global coordinates
3. Performs global frame analysis (stiffness method)
4. Computes nodal displacements and reactions
5. Extracts member end forces

### Properties

```python
@property
def nodal_displacements(self) -> Dict[str, np.ndarray]:
    """
    Dictionary mapping node ID to displacement vector [UX, UY, UZ, RX, RY, RZ].
    Displacements are in global coordinates.
    """

@property
def reactions(self) -> Dict[str, np.ndarray]:
    """
    Dictionary mapping supported node ID to reaction vector [FX, FY, FZ, MX, MY, MZ].
    Reactions are in global coordinates.
    """

@property
def member_end_forces(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Dictionary mapping member ID to (start_forces, end_forces).
    Forces are in local member coordinates [Nx, Vy, Vz, Tx, My, Mz].
    """
```

### Analysis Methods

```python
def get_member_results(self, member_id: str) -> MemberResults:
    """
    Get detailed analysis results for a specific member.
    
    Returns a MemberResults object containing:
    - Axial force distribution
    - Shear force distributions (y and z)
    - Bending moment distributions (y and z)
    - Torsion distribution
    - Displacement distributions
    - Stress distributions
    """

def to_loaded_beams(self) -> Dict[str, LoadedBeam]:
    """
    Convert frame analysis to individual LoadedBeam objects.
    
    Each member is converted to a LoadedBeam with:
    - Member end forces as support reactions
    - Any intermediate member loads applied
    - Local coordinate system orientation
    
    Returns dictionary mapping member ID to LoadedBeam.
    """
```

### MemberResults

Detailed results for an individual member.

```python
@dataclass
class MemberResults:
    member: Member
    
    # Analysis results (all are Result objects like LoadedBeam)
    axial: AnalysisResult          # N(x), sigma_axial(x), u(x)
    shear_y: AnalysisResult        # Vy(x), tau_y(x), v(x) 
    shear_z: AnalysisResult        # Vz(x), tau_z(x), w(x)
    bending_y: AnalysisResult      # My(x), sigma_y(x), theta_y(x)
    bending_z: AnalysisResult      # Mz(x), sigma_z(x), theta_z(x)
    torsion: AnalysisResult        # T(x), tau_torsion(x), phi(x)
    
    @property
    def von_mises(self) -> Result:
        """Combined Von Mises stress distribution."""
```

### Visualization

All 3D plots use matplotlib's wireframe style for clarity. Members are rendered as lines connecting nodes, with optional color mapping for results visualization.

#### Basic Frame Plot

```python
def plot(
    self,
    show_loads: bool = True,
    show_reactions: bool = True,
    show_member_ids: bool = True,
    show_node_ids: bool = True,
    deformed: bool = False,
    scale_factor: float = 1.0,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the frame geometry in 3D wireframe style.
    
    Args:
        show_loads: Display applied load arrows
        show_reactions: Display reaction arrows at supports
        show_member_ids: Label members with their IDs
        show_node_ids: Label nodes with their IDs
        deformed: If True, show deformed shape overlay
        scale_factor: Scale factor for deformed shape visualization
        save_path: Path to save the plot (supports .svg)
    
    Rendering:
        - Undeformed shape: solid black lines
        - Deformed shape: dashed blue lines (when deformed=True)
        - Supports: triangle markers at constrained nodes
        - Loads: arrow quivers in red
        - Reactions: arrow quivers in green
    """
```

#### Deflection Plot

```python
def plot_deflection(
    self,
    scale_factor: float = 1.0,
    points_per_member: int = 20,
    colormap: str = "viridis",
    show_undeformed: bool = True,
    show_colorbar: bool = True,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the deformed frame shape in 3D wireframe style, colored by displacement magnitude.
    
    Args:
        scale_factor: Multiplier for displacement visualization (1.0 = true scale)
        points_per_member: Number of interpolation points along each member
        colormap: Matplotlib colormap name for displacement magnitude
        show_undeformed: If True, show original geometry as faint dashed lines
        show_colorbar: If True, display colorbar with displacement units
        save_path: Path to save the plot (supports .svg)
    
    Rendering:
        - Each member is discretized into segments
        - Segment colors represent total displacement magnitude at that point
        - Displacement magnitude = sqrt(Ux² + Uy² + Uz²)
        - Undeformed shape shown as light gray dashed wireframe
        - Colorbar shows displacement range (min to max)
    
    Example:
        >>> loaded_frame.plot_deflection(
        ...     scale_factor=100,      # Exaggerate deflections 100x
        ...     colormap="plasma",
        ...     save_path="deflection.svg"
        ... )
    """
```

#### Von Mises Stress Plot

```python
def plot_von_mises(
    self,
    points_per_member: int = 20,
    colormap: str = "jet",
    show_colorbar: bool = True,
    stress_limits: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the frame in 3D wireframe style, colored by Von Mises stress.
    
    Args:
        points_per_member: Number of interpolation points along each member
        colormap: Matplotlib colormap name for stress visualization
        show_colorbar: If True, display colorbar with stress units
        stress_limits: Optional (min, max) tuple to fix colorbar range.
                       If None, auto-scales to data range.
        save_path: Path to save the plot (supports .svg)
    
    Rendering:
        - Each member is discretized into line segments
        - Segment colors represent Von Mises stress at that point
        - Von Mises stress computed as: σ_vm = sqrt(σ² + 3τ²)
        - High stress regions shown in warm colors (red/yellow)
        - Low stress regions shown in cool colors (blue/green)
    
    Stress Calculation:
        At each point along each member, Von Mises stress is computed from:
        - σ_axial: Axial stress from N/A
        - σ_bending: Combined bending stress from My and Mz
        - τ_shear: Combined shear stress from Vy and Vz
        - τ_torsion: Torsional shear stress from T
        
        σ_max = |σ_axial| + |σ_bending_y| + |σ_bending_z|  (conservative superposition)
        τ_max = |τ_shear_y| + |τ_shear_z| + |τ_torsion|
        σ_vm = sqrt(σ_max² + 3·τ_max²)
    
    Example:
        >>> loaded_frame.plot_von_mises(
        ...     colormap="hot",
        ...     stress_limits=(0, 250e6),  # 0 to 250 MPa
        ...     save_path="stress.svg"
        ... )
    """
```

#### Combined Stress and Deflection Plot

```python
def plot_results(
    self,
    result_type: str = "von_mises",
    deformed: bool = True,
    scale_factor: float = 1.0,
    points_per_member: int = 20,
    colormap: str = "jet",
    show_undeformed: bool = True,
    show_colorbar: bool = True,
    show_node_ids: bool = False,
    show_member_ids: bool = False,
    value_limits: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> None:
    """
    Unified 3D wireframe plot showing deformed shape colored by analysis results.
    
    Args:
        result_type: Type of result to display. Options:
            - "von_mises": Von Mises stress
            - "deflection": Displacement magnitude
            - "axial_stress": Axial stress (σ from N)
            - "bending_stress": Combined bending stress (σ from M)
            - "shear_stress": Combined shear stress (τ from V)
            - "axial_force": Axial force N
            - "shear_force": Combined shear force magnitude
            - "bending_moment": Combined bending moment magnitude
            - "torsion": Torsional moment T
        deformed: If True, plot on deformed geometry
        scale_factor: Displacement scale factor (only used if deformed=True)
        points_per_member: Number of interpolation points per member
        colormap: Matplotlib colormap name
        show_undeformed: Show original geometry as reference
        show_colorbar: Display colorbar legend
        show_node_ids: Label nodes
        show_member_ids: Label members
        value_limits: Optional (min, max) to fix color range
        save_path: Path to save the plot (supports .svg)
    
    Example:
        >>> # Von Mises stress on deformed shape
        >>> loaded_frame.plot_results(
        ...     result_type="von_mises",
        ...     deformed=True,
        ...     scale_factor=50,
        ...     save_path="results.svg"
        ... )
        
        >>> # Bending moment diagram in 3D
        >>> loaded_frame.plot_results(
        ...     result_type="bending_moment",
        ...     deformed=False,
        ...     colormap="RdBu_r",  # Diverging colormap for +/- values
        ...     save_path="moment.svg"
        ... )
    """
```

#### Member Force Diagrams (2D)

```python
def plot_member_diagrams(
    self,
    member_id: str,
    save_path: Optional[str] = None
) -> None:
    """
    Plot internal force diagrams (N, V, M, T) for a specific member.
    
    Creates a 2x2 subplot figure showing:
        - Top-left: Axial force N(x)
        - Top-right: Shear forces Vy(x) and Vz(x)
        - Bottom-left: Bending moments My(x) and Mz(x)
        - Bottom-right: Torsion T(x)
    
    All diagrams use the member's local x-axis (0 to L).
    """
```

#### Visualization Examples

```python
# Basic geometry plot
loaded_frame.plot(show_loads=True, show_reactions=True)

# Deformed shape with displacement coloring
loaded_frame.plot_deflection(
    scale_factor=100,
    colormap="viridis",
    save_path="deflection.svg"
)

# Von Mises stress distribution
loaded_frame.plot_von_mises(
    colormap="jet",
    stress_limits=(0, 345e6),  # Cap at yield stress
    save_path="von_mises.svg"
)

# Combined: stress on deformed shape
loaded_frame.plot_results(
    result_type="von_mises",
    deformed=True,
    scale_factor=50,
    show_undeformed=True,
    save_path="combined.svg"
)

# Axial force distribution (useful for trusses)
loaded_frame.plot_results(
    result_type="axial_force",
    deformed=False,
    colormap="RdBu_r",  # Red = tension, Blue = compression
    save_path="axial.svg"
)
```

#### Wireframe Rendering Details

The 3D wireframe visualization uses matplotlib's `Axes3D` with `Line3DCollection`:

- **Line Segments**: Each member is divided into `points_per_member` segments
- **Color Mapping**: Segment colors interpolated from result values at endpoints
- **Deformation**: Node positions offset by scaled displacement vectors
- **Interpolation**: Results sampled along member using FEM shape functions

```
Wireframe Structure:
    
    Node B ──────────────── Node C
       │   ╲               ╱   │
       │    ╲   segments  ╱    │
       │     ╲───────────╱     │
       │      colored by       │
       │      result value     │
       │                       │
    Node A                   Node D
    (support)               (support)
```

The wireframe style provides:
- Clear visibility of frame geometry
- Easy interpretation of stress/deflection gradients
- Lightweight rendering for large structures
- Consistent appearance across all plot types

### Example

```python
from beamy.frame import LoadedFrame

# Analyze frame
loaded_frame = LoadedFrame(frame, loads)

# Get reactions at supports
print("Reactions at node A:", loaded_frame.reactions["A"])
print("Reactions at node B:", loaded_frame.reactions["B"])

# Get detailed results for a member
beam_results = loaded_frame.get_member_results("beam1")
print(f"Max bending moment: {beam_results.bending_z.action.abs_max:.2f} N·m")
print(f"Max deflection: {beam_results.shear_z.displacement.abs_max:.6f} m")

# Plot frame
loaded_frame.plot(deformed=True, scale_factor=100)

# Access individual LoadedBeam objects for compatibility
loaded_beams = loaded_frame.to_loaded_beams()
beam1_lb = loaded_beams["beam1"]
beam1_lb.check_aisc_chapter_f("m", "N")  # Run code checks
```

---

## Workflow Example

Complete example showing the frame analysis workflow:

```python
import numpy as np
from sectiony.library import i as i_section
from beamy import Material
from beamy.frame import Node, Member, Frame, FrameLoadCase, LoadedFrame

# ============================================
# 1. DEFINE MATERIALS AND SECTIONS
# ============================================
steel = Material(name="A992 Steel", E=200e9, G=80e9, Fy=345e6)

column_section = i_section(d=0.31, b=0.31, tf=0.017, tw=0.011, r=0.0)  # W310x97
beam_section = i_section(d=0.53, b=0.21, tf=0.016, tw=0.010, r=0.0)    # W530x66

# ============================================
# 2. DEFINE FRAME GEOMETRY  
# ============================================
# Simple portal frame in XZ plane
nodes = [
    Node("A", np.array([0.0, 0.0, 0.0]), support="111111"),   # Fixed base left
    Node("B", np.array([8.0, 0.0, 0.0]), support="111111"),   # Fixed base right
    Node("C", np.array([0.0, 0.0, 5.0])),                     # Top left (free)
    Node("D", np.array([8.0, 0.0, 5.0])),                     # Top right (free)
]

members = [
    # Left column (vertical, y-axis pointing in +X)
    Member("col_left", "A", "C", column_section, steel, np.array([1.0, 0.0, 0.0])),
    # Right column
    Member("col_right", "B", "D", column_section, steel, np.array([1.0, 0.0, 0.0])),
    # Beam (horizontal, y-axis pointing up in +Z for strong-axis bending)
    Member("beam", "C", "D", beam_section, steel, np.array([0.0, 0.0, 1.0])),
]

frame = Frame.from_nodes_and_members(nodes, members)

# ============================================
# 3. DEFINE LOADS
# ============================================
loads = FrameLoadCase("Gravity + Lateral")

# Uniform gravity load on beam (local -z = global -Z for this orientation)
loads.add_member_uniform_force("beam", np.array([0.0, 0.0, -25000.0]))  # 25 kN/m

# Lateral point load at top left (wind)
loads.add_nodal_force("C", np.array([30000.0, 0.0, 0.0]))  # 30 kN in +X

# ============================================
# 4. ANALYZE
# ============================================
loaded_frame = LoadedFrame(frame, loads)

# ============================================
# 5. EXTRACT RESULTS
# ============================================
# Global reactions
print("=== Support Reactions ===")
for node_id, reaction in loaded_frame.reactions.items():
    print(f"Node {node_id}: F={reaction[:3]}, M={reaction[3:]}")

# Member forces
print("\n=== Beam Results ===")
beam = loaded_frame.get_member_results("beam")
print(f"  Max Moment (strong axis): {beam.bending_z.action.abs_max/1000:.1f} kN·m")
print(f"  Max Shear: {beam.shear_z.action.abs_max/1000:.1f} kN")
print(f"  Max Deflection: {beam.shear_z.displacement.abs_max*1000:.2f} mm")

print("\n=== Column Results ===")
col = loaded_frame.get_member_results("col_left")
print(f"  Max Axial: {col.axial.action.abs_max/1000:.1f} kN")
print(f"  Max Moment: {col.bending_z.action.abs_max/1000:.1f} kN·m")

# ============================================
# 6. VISUALIZE (3D Wireframe Plots)
# ============================================
# Basic geometry with loads and reactions
loaded_frame.plot(deformed=True, scale_factor=50, save_path="portal_frame.svg")

# Deflection plot - colored by displacement magnitude
loaded_frame.plot_deflection(
    scale_factor=50,
    colormap="viridis",
    show_undeformed=True,
    save_path="portal_deflection.svg"
)

# Von Mises stress plot - identify high-stress regions
loaded_frame.plot_von_mises(
    colormap="jet",
    stress_limits=(0, steel.Fy),  # Cap colorbar at yield stress
    save_path="portal_von_mises.svg"
)

# Combined: Von Mises stress on deformed shape
loaded_frame.plot_results(
    result_type="von_mises",
    deformed=True,
    scale_factor=50,
    show_undeformed=True,
    colormap="hot",
    save_path="portal_combined.svg"
)

# ============================================
# 7. CODE CHECKS (via LoadedBeam conversion)
# ============================================
loaded_beams = loaded_frame.to_loaded_beams()
for member_id, lb in loaded_beams.items():
    result = lb.check_aisc_chapter_f("m", "N")
    print(f"\n{member_id} AISC Check: {'PASS' if result['pass'] else 'FAIL'}")
```

---

## Design Decisions

### Why Separate Node and Member Definitions?

Explicit node definitions allow:
- Clear support specification at connection points
- Unambiguous node identification for load application
- Easy visualization and debugging
- Support for nodes with multiple connected members (frame joints)

### Why Member-Based Orientation Instead of Rotation Angle?

Using an orientation vector (local y-axis direction) is more intuitive for 3D frames than a rotation angle because:
- No ambiguity about rotation axis or direction
- Works naturally for members in any orientation
- Matches common engineering practice (specifying "web faces this direction")

### Coordinate Transformation Strategy

Loads and results can be expressed in either local or global coordinates:
- **Input flexibility**: Users can specify loads in whichever system is natural
- **Internal consistency**: All analysis done in global coordinates
- **Output flexibility**: Results available in both systems

### Conversion to LoadedBeam Objects

The `to_loaded_beams()` method enables:
- Reuse of existing 1D analysis and plotting infrastructure
- Access to all `LoadedBeam` methods (shear, bending, stress, etc.)
- Compatibility with existing code checks (AISC Chapter F)
- Familiar API for users already using `Beam1D`

### Future Extensions

Planned additions (not in initial implementation):
- Rigid end offsets for eccentric connections
- Temperature loads
- Prestress / initial member forces
- P-delta effects (geometric nonlinearity)
- Member self-weight as automatic load
- Load combinations and envelopes
