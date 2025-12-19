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
    element_type: Literal["beam", "truss", "cable"] = "beam"  # Element behavior type
    releases: Optional[str] = None  # End releases (12-digit string)
    constraints: Optional[str] = None  # Member-end constraints (12-digit string)
```

### Parameters

- `id` (str): Unique member identifier (e.g., "M1", "column_1", "beam_AB")
- `start_node_id` (str): ID of the node at the start of the member
- `end_node_id` (str): ID of the node at the end of the member
- `section` (Section): Cross-section properties (from sectiony library)
- `material` (Material): Material properties (E, G, Fy)
- `orientation` (np.ndarray): A 3D vector in global coordinates that defines the direction of the local y-axis. This vector is projected onto the plane perpendicular to the member axis to determine the actual local y-direction.
- `element_type` (str): Element behavior type (default: "beam")
  - `"beam"`: Full 3D beam element (all DOFs active)
  - `"truss"`: Two-force member (only axial forces, moments released)
  - `"cable"`: Cable element (tension-only, nonlinear)
- `releases` (str, optional): Member end releases as a 12-digit string (see [End Releases](#end-releases))
- `constraints` (str, optional): Member-end rigid constraints as a 12-digit string (see [Member End Constraints](#member-end-constraints))

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
- `"000001000001"`: Pinned at both ends about z-axis (simple beam, free to rotate about z)
- `"000001000000"`: Pinned at start only
- `"111111111111"`: All releases (isolated member - typically not physical)

### Member End Constraints

The `constraints` parameter is a 12-digit string constraining forces/moments at member ends relative to the global frame:
- Digits 1-6: Start node constraints [Nx, Vy, Vz, Tx, My, Mz]
- Digits 7-12: End node constraints [Nx, Vy, Vz, Tx, My, Mz]

Where `0` = free and `1` = constrained (acts like a built-in support).

**Use Case:** Member-end constraints are used for eccentric or special connection conditions that differ from typical frame joints. Not commonly used in standard frame analysis.

### Element Types

#### Beam Element (Default)

Full 3D beam element with all 6 DOFs active at each node. Transmits:
- All internal forces (Nx, Vy, Vz)
- All moments (Tx, My, Mz)

Affected by:
- Bending stiffness (EI)
- Torsional stiffness (GJ)
- Axial stiffness (EA)
- Shear deformation (optional)

```python
Member("beam1", "A", "B", section, steel, orientation, element_type="beam")
```

#### Truss Element

Two-force member that can only transmit axial force. All moments and transverse forces are automatically released.

Characteristics:
- Forces internal to member align with member axis
- Acts like two-force members in trusses
- Typically much smaller rotation stiffness
- Useful for bracing, diagonals

```python
Member("brace1", "A", "C", small_section, steel, orientation, element_type="truss")
```

#### Cable Element

Tension-only element for cable/guy-wire analysis. Cannot transmit compression or bending.

Characteristics:
- Nonlinear behavior (only active in tension)
- Forces must be along member axis (axial only)
- Deactivates when subjected to compression
- Useful for cable-stayed structures, suspension systems

```python
Member("cable1", "A", "D", cable_section, steel, orientation, element_type="cable")
```

### Orientation Examples

```python
# Vertical column - local y points in global +X direction (weak axis perpendicular)
column = Member(
    id="col1",
    start_node_id="A", 
    end_node_id="B",
    section=w_section,
    material=steel,
    orientation=np.array([1.0, 0.0, 0.0]),  # y-axis → global X
    element_type="beam"
)

# Horizontal beam in XZ plane - local y points up (global +Z) for strong-axis bending
beam = Member(
    id="beam1",
    start_node_id="B",
    end_node_id="C", 
    section=w_section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0]),  # y-axis → global Z
    element_type="beam"
)

# Inclined member (brace) - y-axis perpendicular to member, pointing "up"
brace = Member(
    id="brace1",
    start_node_id="A",
    end_node_id="C",
    section=angle_section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0]),  # Projected onto perpendicular plane
    element_type="truss"  # Two-force member (brace)
)

# Cable with releases at both ends (only transmits tension)
cable = Member(
    id="guy_wire",
    start_node_id="D",
    end_node_id="E",
    section=cable_section,
    material=steel,
    orientation=np.array([0.0, 0.0, 1.0]),
    element_type="cable",
    releases="111111111111"  # All forces/moments released (tension only)
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
    member_point_moments: List[MemberPointMoment] = field(default_factory=list)
    member_distributed_forces: List[MemberDistributedForce] = field(default_factory=list)
    nodal_springs: List[NodalSpring] = field(default_factory=list)
```

See the [Loads reference documentation](../setup/loads.md#frameloadcase) for complete `FrameLoadCase` documentation including all load types and methods.

### Quick Reference: Adding Loads

```python
from beamy.frame import FrameLoadCase
import numpy as np

loads = FrameLoadCase(name="Combined Loads")

# Nodal forces (at specific nodes)
loads.add_nodal_force("D", force=np.array([50000.0, 0.0, 0.0]))

# Nodal moments
loads.add_nodal_moment("C", moment=np.array([0.0, 0.0, 25000.0]))

# Point forces along members
loads.add_member_point_force("beam1", position=3.0, force=np.array([0.0, 0.0, -50000.0]))

# Distributed forces over member segments
loads.add_member_distributed_force(
    member_id="beam1",
    start_position=0.0,
    end_position=10.0,
    start_force=np.array([0.0, 0.0, -5000.0]),
    end_force=np.array([0.0, 0.0, -5000.0])
)

# Uniform distributed forces (entire member)
loads.add_member_uniform_force("beam1", force=np.array([0.0, 0.0, -3000.0]))

# Elastic nodal springs (partial boundary conditions)
K = np.diag([1e6, 1e6, 1e6, 1e5, 1e5, 1e5])
loads.add_nodal_spring("A", K)
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
