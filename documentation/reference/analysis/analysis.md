# Analysis

Reference documentation for beam analysis functions and result classes.

## Table of Contents

- [Result](#result)
- [AnalysisResult](#analysisresult)
- [LoadedBeam](#loadedbeam)
- [LoadedFrame](#loadedframe)
- [MemberResults](#memberresults)
- [solve_x_reactions](#solve_x_reactions)
- [solve_transverse_reactions](#solve_transverse_reactions)
- [get_all_loads](#get_all_loads)
- [plot_analysis_results](#plot_analysis_results)

---

## Result

Wrapper class for analysis results providing convenient access to (x, y) data pairs.

```python
from beamy import Result

class Result:
    def __init__(self, x: np.ndarray, values: np.ndarray):
        """
        Args:
            x: Array of x-coordinates
            values: Array of corresponding y-values
        """
```

### Properties

- `max` (float): Maximum value in the result
- `min` (float): Minimum value in the result
- `mean` (float): Mean value
- `range` (float): Range (max - min) of values

### Methods

#### `at(x_loc: float) -> float`

Interpolate the result at a specific x-location.

**Parameters:**
- `x_loc` (float): X-coordinate to evaluate at

**Returns:**
- `float`: Interpolated value at `x_loc`

### Usage

The `Result` class is iterable and can be indexed:

```python
# Iterate over (x, y) pairs
for x, y in result:
    print(f"x={x}, y={y}")

# Access by index
x, y = result[0]

# Get value at specific location
value = result.at(2.5)

# Get statistics
print(f"Max: {result.max}, Min: {result.min}, Mean: {result.mean}")
```

### Example

```python
import numpy as np
from beamy import Result

x = np.linspace(0, 5, 100)
y = np.sin(x)
result = Result(x, y)

print(result.max)      # Maximum value
print(result.at(2.5))  # Interpolated value at x=2.5
```

---

## AnalysisResult

Container for action, stress, and displacement results.

```python
from beamy import AnalysisResult

@dataclass
class AnalysisResult:
    _action: Result
    _stress: Result
    _displacement: Result
```

### Properties

- `action` (Result): Force or moment distribution (V, M, N, T)
- `stress` (Result): Stress distribution (sigma, tau)
- `displacement` (Result): Displacement distribution (w, u, theta)

### Example

```python
# Get analysis result from LoadedBeam
result = loaded_beam.shear(axis="z", points=100)

# Access action (shear force)
shear_force = result.action
print(shear_force.max)  # Maximum shear force

# Access stress (shear stress)
shear_stress = result.stress
print(shear_stress.at(2.5))  # Shear stress at x=2.5

# Access displacement
displacement = result.displacement
print(displacement.min)  # Minimum displacement
```

---

## LoadedBeam

Main analysis class that combines a beam with loads and provides analysis methods.

```python
from beamy import LoadedBeam

@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase
```

### Initialization

When a `LoadedBeam` is created, it automatically:
1. Solves for support reactions in all directions (x, y, z, rotations)
2. Combines applied loads with support reactions
3. Stores the complete load set for analysis

### Methods

#### `shear(axis: str, points: int = 100) -> AnalysisResult`

Calculate shear force and shear stress distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Shear force Fy and bending moment Mz
  - `"z"`: Shear force Fz and bending moment My
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (shear force), stress (shear stress), and displacement

#### `bending(axis: str, points: int = 100) -> AnalysisResult`

Calculate bending moment and bending stress distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Bending in x-z plane (moment My)
  - `"z"`: Bending in x-y plane (moment Mz)
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (bending moment), stress (bending stress), and displacement

#### `axial(points: int = 100) -> AnalysisResult`

Calculate axial force and axial stress distribution.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (axial force N), stress (axial stress σ = N/A), and displacement (axial displacement u)

#### `torsion(points: int = 100) -> AnalysisResult`

Calculate torsional moment and torsional stress distribution.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `AnalysisResult`: Contains action (torsional moment T), stress (torsional stress τ = T·r/J), and displacement (twist angle θ)

#### `deflection(axis: str, points: int = 100) -> Result`

Calculate transverse deflection distribution.

**Parameters:**
- `axis` (str): Bending axis, either `"y"` or `"z"`
  - `"y"`: Deflection in y-direction (v)
  - `"z"`: Deflection in z-direction (w)
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `Result`: Displacement distribution

#### `von_mises(points: int = 100) -> Result`

Calculate the Von Mises stress distribution along the beam.

**Parameters:**
- `points` (int): Number of evaluation points along the beam (default: 100)

**Returns:**
- `Result`: Von Mises stress distribution

**Note:** Uses a conservative superposition of maximum stress components:
- `sigma_vm = sqrt(sigma_max^2 + 3 * tau_max^2)`
- `sigma_max = |sigma_axial| + |sigma_bending_y| + |sigma_bending_z|`
- `tau_max = |tau_shear_y| + |tau_shear_z| + |tau_torsion|`

#### `plot(plot_stress: bool = False, plot_section: bool = True) -> None`

Plot the 3D beam diagram with loads and optional stress visualization.

**Parameters:**
- `plot_stress` (bool): If True, color the beam axis by von Mises stress (default: False)
- `plot_section` (bool): If True, draw the section outline at x=0 (default: True)

**Note:** This is a convenience method that calls `plot_beam_diagram()`. See the [Plotter documentation](plotter.md) for details.

### Example

```python
from beamy import Beam1D, Material, Section, Support, LoadCase, PointForce, LoadedBeam
import numpy as np

# Define beam
steel = Material(name="Steel", E=200e9, G=80e9)
section = Section(
    name="Rectangular",
    A=0.02, Iy=6.67e-5, Iz=1.67e-5, J=1e-6,
    y_max=0.05, z_max=0.1
)
beam = Beam1D(
    L=5.0,
    material=steel,
    section=section,
    supports=[
        Support(x=0.0, type="111000"),
        Support(x=5.0, type="111000")
    ]
)

# Define loads
lc = LoadCase(name="Point Load")
lc.add_point_force(PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10_000.0])
))

# Create loaded beam and analyze
loaded_beam = LoadedBeam(beam, lc)

# Shear analysis
shear_result = loaded_beam.shear(axis="z", points=100)
print(f"Max shear: {shear_result.action.max}")  # force units
print(f"Max shear stress: {shear_result.stress.max}")  # force/length² units

# Bending analysis
bending_result = loaded_beam.bending(axis="z", points=100)
print(f"Max moment: {bending_result.action.max}")  # force×length units
print(f"Max bending stress: {bending_result.stress.max}")  # force/length² units

# Deflection
deflection = loaded_beam.deflection(axis="z", points=100)
print(f"Max deflection: {deflection.max}")  # length units
print(f"Deflection at midspan: {deflection.at(2.5)}")  # length units

# Axial analysis
axial_result = loaded_beam.axial(points=100)

# Torsion analysis
torsion_result = loaded_beam.torsion(points=100)

# Von Mises stress
von_mises = loaded_beam.von_mises(points=100)
print(f"Max von Mises stress: {von_mises.max}")  # force/length² units

# Plot beam diagram
loaded_beam.plot()
```

---

## LoadedFrame

A 3D frame structure with loads applied, ready for analysis.

```python
from beamy.frame import LoadedFrame

@dataclass
class LoadedFrame:
    frame: Frame
    loads: FrameLoadCase
    settings: Optional[FrameAnalysisSettings] = None
```

### Initialization

When a `LoadedFrame` is created, it automatically:
1. Validates load references (node IDs and member IDs exist)
2. Assembles the global stiffness matrix from member properties
3. Applies loads and boundary conditions
4. Solves for nodal displacements and reactions using the direct stiffness method
5. Extracts member end forces and internal stress distributions

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

#### `get_member_results(member_id: str) -> MemberResults`

Get detailed analysis results for a specific member.

**Parameters:**
- `member_id` (str): ID of the member to analyze

**Returns:**
- `MemberResults`: Object containing stress, force, and displacement distributions

#### `to_loaded_beams(self) -> Dict[str, LoadedBeam]`

Convert frame analysis to individual LoadedBeam objects.

Each member is converted to a LoadedBeam with:
- Member end forces as support reactions
- Any intermediate member loads applied
- Local coordinate system orientation

**Returns:**
- `Dict[str, LoadedBeam]`: Dictionary mapping member ID to LoadedBeam object

### Visualization Methods

#### `plot(show_loads: bool = True, show_reactions: bool = True, show_member_ids: bool = True, show_node_ids: bool = True, deformed: bool = False, scale_factor: float = 1.0, save_path: Optional[str] = None) -> None`

Plot the frame geometry in 3D wireframe style.

**Parameters:**
- `show_loads` (bool): Display applied load arrows (default: True)
- `show_reactions` (bool): Display reaction arrows at supports (default: True)
- `show_member_ids` (bool): Label members with their IDs (default: True)
- `show_node_ids` (bool): Label nodes with their IDs (default: True)
- `deformed` (bool): If True, show deformed shape overlay (default: False)
- `scale_factor` (float): Scale factor for deformed shape visualization (default: 1.0)
- `save_path` (str, optional): Path to save the plot (supports .svg)

**Rendering:**
- Undeformed shape: solid black lines
- Deformed shape: dashed blue lines (when deformed=True)
- Supports: triangle markers at constrained nodes
- Loads: arrow quivers in red
- Reactions: arrow quivers in green

#### `plot_deflection(scale_factor: float = 1.0, points_per_member: int = 20, colormap: str = "viridis", show_undeformed: bool = True, show_colorbar: bool = True, save_path: Optional[str] = None) -> None`

Plot the deformed frame shape in 3D wireframe style, colored by displacement magnitude.

**Parameters:**
- `scale_factor` (float): Multiplier for displacement visualization (default: 1.0)
- `points_per_member` (int): Number of interpolation points along each member (default: 20)
- `colormap` (str): Matplotlib colormap name (default: "viridis")
- `show_undeformed` (bool): Show original geometry as faint dashed lines (default: True)
- `show_colorbar` (bool): Display colorbar with displacement range (default: True)
- `save_path` (str, optional): Path to save the plot (supports .svg)

**Rendering:**
- Each member discretized into segments
- Segment colors represent total displacement magnitude (sqrt(Ux² + Uy² + Uz²))
- Undeformed shape shown as light gray dashed wireframe
- Colorbar shows displacement range (min to max)

#### `plot_von_mises(points_per_member: int = 20, colormap: str = "jet", show_colorbar: bool = True, stress_limits: Optional[Tuple[float, float]] = None, save_path: Optional[str] = None) -> None`

Plot the frame in 3D wireframe style, colored by Von Mises stress.

**Parameters:**
- `points_per_member` (int): Number of interpolation points along each member (default: 20)
- `colormap` (str): Matplotlib colormap name (default: "jet")
- `show_colorbar` (bool): Display colorbar with stress range (default: True)
- `stress_limits` (Tuple[float, float], optional): Optional (min, max) tuple to fix colorbar range
- `save_path` (str, optional): Path to save the plot (supports .svg)

**Rendering:**
- Each member discretized into line segments
- Segment colors represent Von Mises stress at that point
- Von Mises stress: σ_vm = sqrt(σ² + 3τ²)
- High stress regions shown in warm colors (red/yellow)
- Low stress regions shown in cool colors (blue/green)

**Stress Calculation:**
At each point along each member, Von Mises stress is computed from:
- σ_axial: Axial stress (N/A)
- σ_bending: Combined bending stress from My and Mz
- τ_shear: Combined shear stress from Vy and Vz
- τ_torsion: Torsional shear stress from T

Conservative superposition: σ_max = |σ_axial| + |σ_bending|, τ_max = |τ_shear| + |τ_torsion|

#### `plot_results(result_type: str = "von_mises", deformed: bool = True, scale_factor: float = 1.0, points_per_member: int = 20, colormap: str = "jet", show_undeformed: bool = True, show_colorbar: bool = True, show_node_ids: bool = False, show_member_ids: bool = False, value_limits: Optional[Tuple[float, float]] = None, save_path: Optional[str] = None) -> None`

Unified 3D wireframe plot showing deformed shape colored by analysis results.

**Parameters:**
- `result_type` (str): Type of result to display (default: "von_mises")
  - "von_mises": Von Mises stress
  - "deflection": Displacement magnitude
  - "axial_stress": Axial stress (σ from N)
  - "bending_stress": Combined bending stress (σ from M)
  - "shear_stress": Combined shear stress (τ from V)
  - "axial_force": Axial force N
  - "shear_force": Combined shear force magnitude
  - "bending_moment": Combined bending moment magnitude
  - "torsion": Torsional moment T
- `deformed` (bool): If True, plot on deformed geometry (default: True)
- `scale_factor` (float): Displacement scale factor (only used if deformed=True) (default: 1.0)
- `points_per_member` (int): Number of interpolation points per member (default: 20)
- `colormap` (str): Matplotlib colormap name (default: "jet")
- `show_undeformed` (bool): Show original geometry as reference (default: True)
- `show_colorbar` (bool): Display colorbar legend (default: True)
- `show_node_ids` (bool): Label nodes (default: False)
- `show_member_ids` (bool): Label members (default: False)
- `value_limits` (Tuple[float, float], optional): Optional (min, max) to fix color range
- `save_path` (str, optional): Path to save the plot (supports .svg)

#### `plot_member_diagrams(member_id: str, save_path: Optional[str] = None) -> None`

Plot internal force diagrams (N, V, M, T) for a specific member.

**Parameters:**
- `member_id` (str): ID of the member to plot
- `save_path` (str, optional): Path to save the plot (supports .svg)

**Creates a 2x2 subplot figure showing:**
- Top-left: Axial force N(x)
- Top-right: Shear forces Vy(x) and Vz(x)
- Bottom-left: Bending moments My(x) and Mz(x)
- Bottom-right: Torsion T(x)

All diagrams use the member's local x-axis (0 to L).

### Example

```python
from beamy.frame import LoadedFrame, Node, Member, Frame, FrameLoadCase
from beamy import Material
from sectiony.library import i as i_section
import numpy as np

# Create frame (see Frame documentation)
steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
column_section = i_section(d=0.3, b=0.3, tf=0.02, tw=0.015, r=0.0)
beam_section = i_section(d=0.4, b=0.2, tf=0.015, tw=0.010, r=0.0)

nodes = [
    Node("A", np.array([0.0, 0.0, 0.0]), support="111111"),
    Node("B", np.array([8.0, 0.0, 0.0]), support="111111"),
    Node("C", np.array([0.0, 0.0, 5.0])),
    Node("D", np.array([8.0, 0.0, 5.0])),
]

members = [
    Member("col_left", "A", "C", column_section, steel, np.array([1.0, 0.0, 0.0])),
    Member("col_right", "B", "D", column_section, steel, np.array([1.0, 0.0, 0.0])),
    Member("beam", "C", "D", beam_section, steel, np.array([0.0, 0.0, 1.0])),
]

frame = Frame.from_nodes_and_members(nodes, members)

# Create loads
loads = FrameLoadCase("Gravity + Lateral")
loads.add_member_uniform_force("beam", np.array([0.0, 0.0, -25000.0]))
loads.add_nodal_force("C", np.array([30000.0, 0.0, 0.0]))

# Analyze
loaded_frame = LoadedFrame(frame, loads)

# Extract results
print("Reactions at A:", loaded_frame.reactions["A"])
print("Reactions at B:", loaded_frame.reactions["B"])

# Member results
beam = loaded_frame.get_member_results("beam")
print(f"Max moment: {beam.bending_z.action.abs_max:.1f} N·m")
print(f"Max shear: {beam.shear_z.action.abs_max:.1f} N")

# Visualization
loaded_frame.plot(deformed=True, scale_factor=100)
loaded_frame.plot_deflection(scale_factor=100, colormap="viridis", save_path="deflection.svg")
loaded_frame.plot_von_mises(colormap="jet", save_path="stress.svg")
loaded_frame.plot_member_diagrams("beam", save_path="diagrams.svg")

# Convert to LoadedBeam objects for code checks
loaded_beams = loaded_frame.to_loaded_beams()
for member_id, lb in loaded_beams.items():
    result = lb.check_aisc_chapter_f("m", "N")
    print(f"{member_id} AISC Check: {'PASS' if result['pass'] else 'FAIL'}")
```

---

## MemberResults

Detailed analysis results for an individual frame member.

```python
from beamy.frame import MemberResults

@dataclass
class MemberResults:
    member: Member
    
    # Analysis results (all are AnalysisResult objects)
    axial: AnalysisResult              # N(x), sigma_axial(x), u(x)
    shear_y: AnalysisResult            # Vy(x), tau_y(x), v(x)
    shear_z: AnalysisResult            # Vz(x), tau_z(x), w(x)
    bending_y: AnalysisResult          # My(x), sigma_y(x), theta_y(x)
    bending_z: AnalysisResult          # Mz(x), sigma_z(x), theta_z(x)
    torsion: AnalysisResult            # T(x), tau_torsion(x), phi(x)
```

### Properties

```python
@property
def von_mises(self) -> Result:
    """Combined Von Mises stress distribution along the member."""
```

### Analysis Results

Each analysis result (e.g., `axial`, `bending_z`) is an `AnalysisResult` object containing:

- **action**: Force or moment distribution (force/moment units)
  - `axial.action`: Axial force N(x)
  - `shear_z.action`: Shear force Vz(x)
  - `bending_z.action`: Bending moment Mz(x)
  - `torsion.action`: Torsional moment T(x)

- **stress**: Stress distribution (stress units = force/length²)
  - `axial.stress`: Axial stress σ = N/A
  - `shear_z.stress`: Shear stress τ_z
  - `bending_z.stress`: Bending stress σ_z = Mz/Iz
  - `torsion.stress`: Torsional stress τ_t = T·r/J

- **displacement**: Displacement or rotation distribution (displacement/radian units)
  - `axial.displacement`: Axial displacement u(x)
  - `shear_z.displacement`: Transverse displacement w(x)
  - `bending_z.displacement`: Rotation θ_z(x)
  - `torsion.displacement`: Twist angle φ(x)

### Accessing Results

```python
# Get member results from a LoadedFrame
results = loaded_frame.get_member_results("beam1")

# Access action (force/moment) distribution
max_moment = results.bending_z.action.max
moment_at_midspan = results.bending_z.action.at(L/2)

# Access stress distribution
max_stress = results.bending_z.stress.max
stress_at_x = results.bending_z.stress.at(3.0)

# Access displacement distribution
max_deflection = results.shear_z.displacement.max
deflection_history = results.shear_z.displacement

# Access Von Mises stress
vm_max = results.von_mises.max
vm_at_point = results.von_mises.at(2.5)

# Iterate over results
for x, moment_value in results.bending_z.action:
    print(f"At x={x}: M={moment_value}")
```

### Example

```python
loaded_frame = LoadedFrame(frame, loads)
beam = loaded_frame.get_member_results("beam1")

print(f"Max Axial Force: {beam.axial.action.abs_max:.1f} N")
print(f"Max Axial Stress: {beam.axial.stress.abs_max:.2f} Pa")
print(f"Max Axial Displacement: {beam.axial.displacement.abs_max:.6f} m")

print(f"Max Shear (Z): {beam.shear_z.action.abs_max:.1f} N")
print(f"Max Shear Stress (Z): {beam.shear_z.stress.abs_max:.2f} Pa")

print(f"Max Bending Moment (Z): {beam.bending_z.action.abs_max:.1f} N·m")
print(f"Max Bending Stress (Z): {beam.bending_z.stress.abs_max:.2f} Pa")
print(f"Max Deflection (Z): {beam.shear_z.displacement.abs_max:.6f} m")

print(f"Max Torsion: {beam.torsion.action.abs_max:.1f} N·m")
print(f"Max Torsional Stress: {beam.torsion.stress.abs_max:.2f} Pa")

print(f"Max Von Mises Stress: {beam.von_mises.abs_max:.2f} Pa")
```

Solve for axial and torsional reactions using 1D finite element method.

```python
from beamy import solve_x_reactions

def solve_x_reactions(
    beam: Beam1D,
    loads: LoadCase
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve for axial (Fx) and torsional (Mx) reactions.
    
    Args:
        beam: Beam1D object
        loads: LoadCase containing applied loads
        
    Returns:
        Tuple of (d_x, d_rx, x_nodes_axial, x_nodes_torsion)
    """
```

### Parameters

- `beam` (Beam1D): Beam object with supports
- `loads` (LoadCase): Load case containing applied loads

### Returns

- `Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]`: 
  - `d_x`: Global axial displacements
  - `d_rx`: Global torsional displacements
  - `x_nodes_axial`: Node positions for axial FEM
  - `x_nodes_torsion`: Node positions for torsion FEM

### Behavior

- Uses 1D finite element method with linear shape functions
- Automatically updates support reactions in the `Support.reactions` dictionary
- Reactions are stored with keys: `"Fx"` for axial reactions and `"Mx"` for torsional reactions

### Example

```python
from beamy import solve_x_reactions

# Solve for reactions
d_x, d_rx, x_nodes_a, x_nodes_t = solve_x_reactions(beam, load_case)

# Reactions are automatically stored in support.reactions
for support in beam.supports:
    print(f"At x={support.x}: Fx={support.reactions['Fx']}, Mx={support.reactions['Mx']}")
```

---

## solve_transverse_reactions

Solve for transverse reactions (shear and bending) using 1D finite element method with Hermite shape functions.

```python
from beamy import solve_transverse_reactions

def solve_transverse_reactions(
    beam: Beam1D,
    loads: LoadCase,
    axis: str = "z"
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve for transverse reactions (shear and bending moment).
    
    Args:
        beam: Beam1D object
        loads: LoadCase containing applied loads
        axis: Bending axis, either "y" or "z" (default: "z")
        
    Returns:
        Tuple of (d, x_nodes)
    """
```

### Parameters

- `beam` (Beam1D): Beam object with supports
- `loads` (LoadCase): Load case containing applied loads
- `axis` (str): Bending axis, either `"y"` or `"z"` (default: `"z"`)
  - `"y"`: Bending in x-y plane (shear Fy, moment Mz)
  - `"z"`: Bending in x-z plane (shear Fz, moment My)

### Returns

- `Tuple[np.ndarray, np.ndarray]`:
  - `d`: Global displacement vector containing [w1, θ1, w2, θ2, ...]
  - `x_nodes`: Array of node x-positions used in the FEM solution

### Behavior

- Uses 1D finite element method with Hermite cubic shape functions
- Automatically updates support reactions in the `Support.reactions` dictionary
- For `axis="z"`: Reactions stored as `"Fz"` and `"My"`
- For `axis="y"`: Reactions stored as `"Fy"` and `"Mz"`

### Example

```python
from beamy import solve_transverse_reactions

# Solve for reactions in z-direction (vertical)
d_z, x_nodes = solve_transverse_reactions(beam, load_case, axis="z")

# Reactions are automatically stored
for support in beam.supports:
    print(f"At x={support.x}: Fz={support.reactions['Fz']}, My={support.reactions['My']}")

# Solve for reactions in y-direction (horizontal)
d_y, x_nodes = solve_transverse_reactions(beam, load_case, axis="y")
```

---

## get_all_loads

Get a sorted list of all loads (applied + reactions) for analysis.

```python
from beamy import get_all_loads

def get_all_loads(
    loads: LoadCase,
    beam: Beam1D
) -> List[Tuple[float, str, float]]:
    """
    Returns a sorted list of all loads and support reactions.
    
    Args:
        loads: LoadCase containing applied loads
        beam: Beam1D object with supports (reactions should be computed)
        
    Returns:
        List of (x, type, magnitude) tuples where:
        - x: Position along beam
        - type: Load type ("Fx", "Fy", "Fz", "Mx", "My", "Mz", "Rx", "Ry", "Rz", "RMx", "RMy", "RMz")
        - magnitude: Load magnitude
    """
```

### Parameters

- `loads` (LoadCase): Load case containing applied loads
- `beam` (Beam1D): Beam object with supports (reactions must be computed first)

### Returns

- `List[Tuple[float, str, float]]`: Sorted list of (x, type, magnitude) tuples

### Load Types

- **Applied Forces**: `"Fx"`, `"Fy"`, `"Fz"`
- **Applied Moments**: `"Mx"`, `"My"`, `"Mz"`
- **Reactions**: `"Rx"`, `"Ry"`, `"Rz"` (force reactions), `"RMx"`, `"RMy"`, `"RMz"` (moment reactions)

### Behavior

- Combines applied loads from the `LoadCase` with support reactions from `beam.supports`
- Loads at the same position and of the same type are summed
- Results are sorted by x-position

### Example

```python
from beamy import get_all_loads, LoadedBeam

# Create loaded beam (reactions are computed automatically)
loaded_beam = LoadedBeam(beam, load_case)

# Get all loads (applied + reactions)
all_loads = get_all_loads(load_case, beam)

# Print all loads
for x, load_type, magnitude in all_loads:
    print(f"At x={x:.2f}: {load_type} = {magnitude:.2f}")

# Filter for specific load type
vertical_forces = [(x, t, v) for x, t, v in all_loads if t in ("Fz", "Rz")]
```

---

## plot_analysis_results

Plots the analysis results (Shear, Moment, Deflection, Axial/Torsion) on 2D line graphs.

```python
from beamy import plot_analysis_results

def plot_analysis_results(
    loaded_beam: LoadedBeam, 
    save_path: Optional[str] = None, 
    show: bool = True, 
    points: int = 100
) -> None:
    """
    Plots the analysis results (Shear, Moment, Deflection, Axial/Torsion) on 2D line graphs.
    """
```

### Parameters

- `loaded_beam` (LoadedBeam): The analyzed LoadedBeam object.
- `save_path` (str, optional): Path to save the figure.
- `show` (bool): Whether to display the plot.
- `points` (int): Number of points to sample.

### Example

```python
plot_analysis_results(loaded_beam)
```

---

## Analysis Workflow

Typical workflow for analyzing a beam:

1. **Define beam geometry and properties**
   ```python
   beam = Beam1D(L=5.0, material=steel, section=section, supports=[...])
   ```

2. **Define loads**
   ```python
   lc = LoadCase(name="My Loads")
   lc.add_point_force(...)
   lc.add_distributed_force(...)
   ```

3. **Create LoadedBeam** (automatically solves for reactions)
   ```python
   loaded_beam = LoadedBeam(beam, lc)
   ```

4. **Perform analysis**
   ```python
   shear_result = loaded_beam.shear(axis="z", points=100)
   bending_result = loaded_beam.bending(axis="z", points=100)
   deflection = loaded_beam.deflection(axis="z", points=100)
   ```

5. **Access results**
   ```python
   max_shear = shear_result.action.max
   max_stress = bending_result.stress.max
   max_deflection = deflection.max
   value_at_x = shear_result.action.at(2.5)
   ```

---

## Units

**Beamy is unit-agnostic.** All results use the same units as your inputs. See the [Beam documentation](../beam.md#units) for details on unit consistency requirements.

- **Action results** (force, moment): Same units as input forces/moments
- **Stress results**: Force / Length² (consistent with your input units)
- **Displacement results**: Length units (same as input length)
- **Rotation results**: Radians (dimensionless)
