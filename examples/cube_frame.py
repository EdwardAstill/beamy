"""
3D Cube Frame Example

This example demonstrates a true 3D frame structure:
- Cube outline (3m × 3m × 3m)
- All 4 base nodes fixed
- Force applied at top corner pointing towards center of cube
- 12 members: 4 vertical columns + 4 base beams + 4 top beams

This shows:
1. True 3D geometry (all three dimensions active)
2. Multiple supports
3. Force with components in all three directions
4. 3D frame behavior with torsion and biaxial bending
"""

import numpy as np
from pathlib import Path
from sectiony.library import i as i_section

from beamy import (
    LoadCase,
    Material,
    NodalForce,
    Frame,
    Member,
    plot_deflection,
    plot_frame,
    plot_von_mises,
)
from beamy.frame import FrameBuilder

# Create gallery subdirectories
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)
frame_dir = gallery_dir / "frame"
frame_dir.mkdir(exist_ok=True)

# ============================================
# 1. MATERIALS AND SECTIONS
# ============================================
steel = Material(name="Steel", E=200e9, G=80e9, Fy=250e6)

# Use square tube sections for the cube members
section = i_section(d=0.2, b=0.2, tf=0.015, tw=0.010, r=0.0)

# ============================================
# 2. DEFINE CUBE GEOMETRY
# ============================================
L = 3.0  # Cube side length (3 meters)
H = 3.0  # Cube height (3 meters)

# Use FrameBuilder for cleaner coordinate-based definition
fb = FrameBuilder(default_orientation=(0, 0, 1))

# Vertical columns (4)
fb.add("col_AE", (0.0, 0.0, 0.0), (0.0, 0.0, H), section, steel, orientation=(1.0, 0.0, 0.0))  # Front-left
fb.add("col_BF", (L, 0.0, 0.0), (L, 0.0, H), section, steel, orientation=(0.0, 1.0, 0.0))      # Front-right
fb.add("col_CG", (L, L, 0.0), (L, L, H), section, steel, orientation=(1.0, 0.0, 0.0))          # Back-right
fb.add("col_DH", (0.0, L, 0.0), (0.0, L, H), section, steel, orientation=(0.0, 1.0, 0.0))      # Back-left

# Base beams (4) - connecting base nodes
fb.add("base_AB", (0.0, 0.0, 0.0), (L, 0.0, 0.0), section, steel)  # Front
fb.add("base_BC", (L, 0.0, 0.0), (L, L, 0.0), section, steel)      # Right
fb.add("base_CD", (L, L, 0.0), (0.0, L, 0.0), section, steel)      # Back
fb.add("base_DA", (0.0, L, 0.0), (0.0, 0.0, 0.0), section, steel)  # Left

# Top beams (4) - connecting top nodes
fb.add("top_EF", (0.0, 0.0, H), (L, 0.0, H), section, steel)       # Front
fb.add("top_FG", (L, 0.0, H), (L, L, H), section, steel)           # Right
fb.add("top_GH", (L, L, H), (0.0, L, H), section, steel)           # Back
fb.add("top_HE", (0.0, L, H), (0.0, 0.0, H), section, steel)       # Left

# Apply fixed supports at all base nodes (z=0)
fb.support_at((0.0, 0.0, 0.0), "111111")
fb.support_at((L, 0.0, 0.0), "111111")
fb.support_at((L, L, 0.0), "111111")
fb.support_at((0.0, L, 0.0), "111111")

# Build frame
frame, coord_to_node = fb.build_with_node_map()

print("=" * 70)
print("3D CUBE FRAME GEOMETRY")
print("=" * 70)
print(f"Cube dimensions: {L}m × {L}m × {H}m")
print(f"Frame: {len(frame.nodes)} nodes, {len(frame.members)} members")
print(f"\nBase nodes (fixed): A, B, C, D")
print(f"Top nodes (free): E, F, G, H")

# ============================================
# 3. APPLY FORCE AT TOP CORNER
# ============================================
loads = LoadCase(name="Force at top corner towards center")

# Force at node G (top corner at [L, L, H])
# Direction: towards center of cube at [L/2, L/2, H/2]
corner_pos = np.array([L, L, H])
center_pos = np.array([L/2, L/2, H/2])
force_direction = center_pos - corner_pos  # [-L/2, -L/2, -H/2]

# Normalize and scale to desired magnitude
force_magnitude = 50000.0  # 50 kN
force_direction_normalized = force_direction / np.linalg.norm(force_direction)
force_vector = force_magnitude * force_direction_normalized

print(f"\n{'=' * 70}")
print("APPLIED LOAD")
print("=" * 70)
print(f"Node: {top_corner_node} (top corner at [{L:.1f}, {L:.1f}, {H:.1f}])")
print(f"Center of cube: [{L/2:.1f}, {L/2:.1f}, {H/2:.1f}]")
print(f"Force magnitude: {force_magnitude/1000:.1f} kN")
print(f"Force direction (normalized): [{force_direction_normalized[0]:.3f}, {force_direction_normalized[1]:.3f}, {force_direction_normalized[2]:.3f}]")
print(f"Force components:")
print(f"  FX = {force_vector[0]/1000:.2f} kN (towards -X)")
print(f"  FY = {force_vector[1]/1000:.2f} kN (towards -Y)")
print(f"  FZ = {force_vector[2]/1000:.2f} kN (downward)")

# Get node ID at top corner (L, L, H)
top_corner_node = coord_to_node[(L, L, H)]
loads.add_nodal_force(top_corner_node, force_vector)

# ============================================
# 4. ANALYZE
# ============================================
print(f"\n{'=' * 70}")
print("PERFORMING 3D FRAME ANALYSIS")
print("=" * 70)

frame.analyze(loads)

print("Analysis complete!")

# ============================================
# 5. RESULTS
# ============================================
print(f"\n{'=' * 70}")
print("SUPPORT REACTIONS")
print("=" * 70)

total_reaction = np.zeros(3)
# Get base node IDs
base_nodes = [coord_to_node[(0.0, 0.0, 0.0)], coord_to_node[(L, 0.0, 0.0)],
              coord_to_node[(L, L, 0.0)], coord_to_node[(0.0, L, 0.0)]]

for node_id in base_nodes:
    if node_id in frame.reactions:
        reaction = frame.reactions[node_id]
    else:
        reaction = np.zeros(6, dtype=float)
    forces = reaction[:3]
    moments = reaction[3:]
    total_reaction += forces
    
    print(f"\nNode {node_id}:")
    print(f"  Forces:  FX={forces[0]/1000:7.2f} kN, FY={forces[1]/1000:7.2f} kN, FZ={forces[2]/1000:7.2f} kN")
    print(f"  Moments: MX={moments[0]/1000:7.2f} kN·m, MY={moments[1]/1000:7.2f} kN·m, MZ={moments[2]/1000:7.2f} kN·m")

print(f"\nTotal reaction forces:")
print(f"  FX = {total_reaction[0]/1000:7.2f} kN (should equal {force_vector[0]/1000:.2f} kN)")
print(f"  FY = {total_reaction[1]/1000:7.2f} kN (should equal {force_vector[1]/1000:.2f} kN)")
print(f"  FZ = {total_reaction[2]/1000:7.2f} kN (should equal {force_vector[2]/1000:.2f} kN)")

# Check equilibrium
print(f"\nEquilibrium check:")
print(f"  X-direction: {abs(total_reaction[0] + force_vector[0]) < 1.0} (error: {abs(total_reaction[0] + force_vector[0]):.2e} N)")
print(f"  Y-direction: {abs(total_reaction[1] + force_vector[1]) < 1.0} (error: {abs(total_reaction[1] + force_vector[1]):.2e} N)")
print(f"  Z-direction: {abs(total_reaction[2] + force_vector[2]) < 1.0} (error: {abs(total_reaction[2] + force_vector[2]):.2e} N)")

# ============================================
# 6. NODAL DISPLACEMENTS
# ============================================
print(f"\n{'=' * 70}")
print("TOP NODE DISPLACEMENTS")
print("=" * 70)

for node_id in ["E", "F", "G", "H"]:
    disp = frame.nodal_displacements[node_id]
    disp_mag = np.linalg.norm(disp[:3])
    
    print(f"\nNode {node_id}:")
    print(f"  Displacement: UX={disp[0]*1000:7.2f} mm, UY={disp[1]*1000:7.2f} mm, UZ={disp[2]*1000:7.2f} mm")
    print(f"  Magnitude: {disp_mag*1000:.2f} mm")
    print(f"  Rotations: RX={disp[3]:7.4f} rad, RY={disp[4]:7.4f} rad, RZ={disp[5]:7.4f} rad")

# ============================================
# 7. MEMBER STRESSES
# ============================================
print(f"\n{'=' * 70}")
print("CRITICAL MEMBER RESULTS")
print("=" * 70)

# Find most stressed members
max_stress = 0
max_stress_member = None

print("\nColumn stresses:")
for col_id in ["col_AE", "col_BF", "col_CG", "col_DH"]:
    profile = frame.demand_provider.actions(col_id, points=401)
    member = frame.get_member(col_id)

    A = member.section.A
    Iy = member.section.Iy
    Iz = member.section.Iz
    J = member.section.J
    y_max = member.section.y_max
    z_max = member.section.z_max
    r_max = max(abs(y_max), abs(z_max))

    N = profile.axial._values
    Vy = profile.shear_y._values
    Vz = profile.shear_z._values
    T = profile.torsion._values
    My = profile.bending_y._values
    Mz = profile.bending_z._values

    sigma_axial = N / A if A > 0 else np.zeros_like(N)
    sigma_bending_y = np.abs(My) * z_max / Iy if Iy > 0 else np.zeros_like(My)
    sigma_bending_z = np.abs(Mz) * y_max / Iz if Iz > 0 else np.zeros_like(Mz)
    sigma = np.abs(sigma_axial) + sigma_bending_y + sigma_bending_z

    tau_shear = (np.abs(Vy) + np.abs(Vz)) / A if A > 0 else np.zeros_like(Vy)
    tau_torsion = np.abs(T) * r_max / J if J > 0 else np.zeros_like(T)
    tau = tau_shear + tau_torsion

    vm = np.sqrt(sigma**2 + 3.0 * tau**2)
    vm_max = float(np.max(vm))
    print(f"  {col_id}: {vm_max/1e6:.1f} MPa")
    
    if vm_max > max_stress:
        max_stress = vm_max
        max_stress_member = col_id

print(f"\nMost stressed member: {max_stress_member} with {max_stress/1e6:.1f} MPa")
print(f"Yield stress: {steel.Fy/1e6:.1f} MPa")
print(f"Utilization: {max_stress/steel.Fy*100:.1f}%")

# ============================================
# 8. VISUALIZATION
# ============================================
print(f"\n{'=' * 70}")
print("GENERATING 3D VISUALIZATIONS")
print("=" * 70)

# Basic geometry with force and reactions
print("1. Cube frame geometry with loads...")
plot_frame(
    frame,
    deformed=True,
    scale_factor=100,
    show_loads=True,
    show_reactions=True,
    save_path=str(frame_dir / "cube_frame_geometry.svg")
)

# Deflection visualization
print("2. Deflection pattern...")
plot_deflection(
    frame,
    scale_factor=100,
    colormap="viridis",
    show_undeformed=True,
    save_path=str(frame_dir / "cube_frame_deflection.svg")
)

# Von Mises stress
print("3. Von Mises stress distribution...")
plot_von_mises(
    frame,
    colormap="turbo",
    stress_limits=(0, steel.Fy),
    save_path=str(frame_dir / "cube_frame_von_mises.svg")
)

print(f"\n{'=' * 70}")
print("COMPLETE!")
print("=" * 70)
print(f"\nPlots saved to: {frame_dir.absolute()}/")
print("  - cube_frame_geometry.svg")
print("  - cube_frame_deflection.svg")
print("  - cube_frame_von_mises.svg")

print(f"\n{'=' * 70}")
print("KEY FEATURES DEMONSTRATED:")
print("=" * 70)
print("* True 3D frame structure (cube outline)")
print("* Multiple fixed supports (all 4 base nodes)")
print("* 3D force vector (components in X, Y, and Z)")
print("* Force direction towards interior point")
print("* 3D structural behavior with torsion and biaxial bending")
print("* Complex load path through spatial frame")
print("* 3D wireframe visualization")

