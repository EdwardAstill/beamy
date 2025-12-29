"""
Portal Frame Example - Demonstrating 3D Frame Analysis

This example shows a 2D portal frame (in the XZ plane) with:
- Fixed supports at the base
- Gravity load on the beam (distributed load applied to member)
- Lateral wind load at top left node (nodal load)
- Point load at beam midspan (member point load)

Demonstrates:
1. Frame geometry definition (nodes and members)
2. Multiple load types:
   - Nodal forces (wind load)
   - Member distributed forces (uniform gravity load)
   - Member point forces (concentrated load at midspan)
3. 3D frame analysis with direct stiffness method
4. Result extraction (reactions, displacements, member forces)
5. Visualization (geometry, deflection, stress)
"""

import numpy as np
from pathlib import Path
from sectiony.library import i as i_section

# Frame-first API
from beamy import (
    LoadCase,
    Material,
    MemberDistributedForce,
    MemberPointForce,
    NodalForce,
    Frame,
    Member,
    Node,
    plot_deflection,
    plot_frame,
    plot_member_diagrams,
    plot_von_mises,
)

# Create gallery subdirectories
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)
frame_dir = gallery_dir / "frame"
frame_dir.mkdir(exist_ok=True)

# ============================================
# 1. DEFINE MATERIALS AND SECTIONS
# ============================================
steel = Material(name="A992 Steel", E=200e9, G=80e9, Fy=345e6)

# W-shape columns and beams
column_section = i_section(d=0.31, b=0.31, tf=0.017, tw=0.011, r=0.0)  # W310x97
beam_section = i_section(d=0.53, b=0.21, tf=0.016, tw=0.010, r=0.0)    # W530x66

# ============================================
# 2. DEFINE FRAME GEOMETRY
# ============================================
# Simple 2D portal frame in XZ plane
# X = horizontal, Z = vertical

members = [
    # Left column (vertical from base to top left)
    # Orientation: local y-axis points in +X direction (strong axis bending for lateral loads)
    Member("col_left", 
           start=np.array([0.0, 0.0, 0.0]),
           end=np.array([0.0, 0.0, 5.0]),
           section=column_section,
           material=steel,
           orientation=np.array([1.0, 0.0, 0.0])),
    
    # Right column (vertical from base to top right)
    Member("col_right",
           start=np.array([8.0, 0.0, 0.0]),
           end=np.array([8.0, 0.0, 5.0]),
           section=column_section,
           material=steel,
           orientation=np.array([1.0, 0.0, 0.0])),
    
    # Beam (horizontal from top left to top right)
    # Orientation: local y-axis points in +Y direction
    # This makes local z point in +Z, so loads in local -z go to global -Z (gravity)
    Member("beam",
           start=np.array([0.0, 0.0, 5.0]),
           end=np.array([8.0, 0.0, 5.0]),
           section=beam_section,
           material=steel,
           orientation=np.array([0.0, 1.0, 0.0])),
]

# Create frame (nodes auto-generated from member endpoints)
frame = Frame.from_members(members)

# Apply supports to the base nodes (auto-generated as N0 and N1)
# Find nodes at base (z=0)
for node in frame.nodes.values():
    if abs(node.position[2]) < 1e-6:  # z ≈ 0
        node.support = "111111"  # Fixed base

print("=" * 60)
print("PORTAL FRAME GEOMETRY")
print("=" * 60)
print(f"Frame: {len(frame.nodes)} nodes, {len(frame.members)} members")
print(f"\nMember lengths:")
for member_id, length in frame.member_lengths.items():
    print(f"  {member_id}: {length:.2f} m")

# ============================================
# 3. DEFINE LOADS - DEMONSTRATING ALL LOAD TYPES
# ============================================
loads = LoadCase(name="Gravity + Lateral + Point Load")

# LOAD TYPE 1: Uniform distributed load on beam (member load)
# This is a gravity load distributed along the entire beam length
# Applied in LOCAL coordinates: local -z direction (downward for this beam)
print(f"\n{'=' * 60}")
print("APPLIED LOADS")
print("=" * 60)

gravity_load = 25000.0  # 25 kN/m
print(f"1. Distributed load on beam: {gravity_load/1000:.1f} kN/m (local -z, gravity)")
loads.add_member_uniform_force("beam", np.array([0.0, 0.0, -gravity_load]))

# LOAD TYPE 2: Point load at beam midspan (member point load)
# Simulating a concentrated load at the center of the beam
midspan_load = 50000.0  # 50 kN
beam_length = frame.get_member("beam").length
print(f"2. Point load at beam midspan: {midspan_load/1000:.1f} kN (local -z)")
loads.add_member_point_force(
    "beam",
    position=0.5,  # 50% along the beam
    force=np.array([0.0, 0.0, -midspan_load]),
    coords="local",
    position_type="relative"  # Position as fraction (0 to 1)
)

# LOAD TYPE 3: Lateral wind load at top left node (nodal load)
# Applied in GLOBAL coordinates: +X direction (wind from left)
# Find top left node (x=0, z=5)
top_left_node_id = None
for nid, node in frame.nodes.items():
    if abs(node.position[0]) < 1e-6 and abs(node.position[2] - 5.0) < 1e-6:
        top_left_node_id = nid
        break

wind_load = 30000.0  # 30 kN
print(f"3. Lateral wind load at node {top_left_node_id}: {wind_load/1000:.1f} kN (global +X)")
loads.add_nodal_force(top_left_node_id, np.array([wind_load, 0.0, 0.0]))

# ============================================
# 4. ANALYZE FRAME
# ============================================
print(f"\n{'=' * 60}")
print("PERFORMING FRAME ANALYSIS")
print("=" * 60)
print("Using direct stiffness method...")

frame.analyze(loads)

print("Analysis complete!")

# ============================================
# 5. EXTRACT RESULTS
# ============================================
print(f"\n{'=' * 60}")
print("SUPPORT REACTIONS")
print("=" * 60)

total_reaction = np.zeros(3)
for node_id, reaction in frame.reactions.items():
    forces = reaction[:3]  # [FX, FY, FZ]
    moments = reaction[3:]  # [MX, MY, MZ]
    total_reaction += forces
    
    print(f"\nNode {node_id}:")
    print(f"  Forces:  FX={forces[0]/1000:8.2f} kN, FY={forces[1]/1000:8.2f} kN, FZ={forces[2]/1000:8.2f} kN")
    print(f"  Moments: MX={moments[0]/1000:8.2f} kN·m, MY={moments[1]/1000:8.2f} kN·m, MZ={moments[2]/1000:8.2f} kN·m")

print(f"\nTotal reaction forces:")
print(f"  FX={total_reaction[0]/1000:8.2f} kN (should = {wind_load/1000:.2f} kN)")
print(f"  FY={total_reaction[1]/1000:8.2f} kN (should = 0.00 kN)")
print(f"  FZ={total_reaction[2]/1000:8.2f} kN (should = {(gravity_load * beam_length + midspan_load)/1000:.2f} kN)")

# Verify equilibrium
applied_vertical = gravity_load * beam_length + midspan_load
applied_horizontal = wind_load
print(f"\nEquilibrium check:")
print(f"  Horizontal: {abs(total_reaction[0] + applied_horizontal) < 1e-6} (error: {abs(total_reaction[0] + applied_horizontal):.2e} N)")
print(f"  Vertical:   {abs(total_reaction[2] - applied_vertical) < 1e-6} (error: {abs(total_reaction[2] - applied_vertical):.2e} N)")

# ============================================
# 6. MEMBER RESULTS
# ============================================
print(f"\n{'=' * 60}")
print("MEMBER ANALYSIS RESULTS")
print("=" * 60)

# Beam results
print("\nBEAM (horizontal member):")
beam_profile = frame.demand_provider.actions("beam", points=801)
print(f"  Max bending moment (about local y): {beam_profile.bending_y.abs_max/1000:.1f} kN·m")
print(f"  Max shear force (local z): {beam_profile.shear_z.abs_max/1000:.1f} kN")

# Column results
print("\nLEFT COLUMN:")
col_left_profile = frame.demand_provider.actions("col_left", points=401)
print(f"  Max axial force: {col_left_profile.axial.abs_max/1000:.1f} kN")
print(f"  Max bending moment (about local y): {col_left_profile.bending_y.abs_max/1000:.1f} kN·m")
print(f"  Max bending moment (about local z): {col_left_profile.bending_z.abs_max/1000:.1f} kN·m")

print("\nRIGHT COLUMN:")
col_right_profile = frame.demand_provider.actions("col_right", points=401)
print(f"  Max axial force: {col_right_profile.axial.abs_max/1000:.1f} kN")
print(f"  Max bending moment (about local y): {col_right_profile.bending_y.abs_max/1000:.1f} kN·m")
print(f"  Max bending moment (about local z): {col_right_profile.bending_z.abs_max/1000:.1f} kN·m")

# ============================================
# 7. VISUALIZE - 3D WIREFRAME PLOTS
# ============================================
print(f"\n{'=' * 60}")
print("GENERATING VISUALIZATION PLOTS")
print("=" * 60)

# Basic geometry with loads and reactions
print("1. Frame geometry with loads and reactions...")
plot_frame(
    frame,
    deformed=True,
    scale_factor=50,
    save_path=str(frame_dir / "portal_frame_geometry.svg")
)

# Deflection plot - colored by displacement magnitude
print("2. Deflection plot (displacement magnitude)...")
plot_deflection(
    frame,
    scale_factor=50,
    colormap="viridis",
    show_undeformed=True,
    save_path=str(frame_dir / "portal_frame_deflection.svg")
)

# Von Mises stress plot
print("3. Von Mises stress distribution...")
plot_von_mises(
    frame,
    colormap="turbo",
    stress_limits=(0, steel.Fy),  # Cap colorbar at yield stress
    save_path=str(frame_dir / "portal_frame_von_mises.svg")
)

# Member force diagrams
print("4. Beam internal force diagrams...")
plot_member_diagrams(
    frame,
    "beam",
    save_path=str(frame_dir / "portal_frame_beam_diagrams.svg")
)

print(f"\n{'=' * 60}")
print("COMPLETE!")
print("=" * 60)
print(f"\nPlots saved to: {frame_dir.absolute()}/")
print("  - portal_frame_geometry.svg")
print("  - portal_frame_deflection.svg")
print("  - portal_frame_von_mises.svg")
print("  - portal_frame_beam_diagrams.svg")

print(f"\n{'=' * 60}")
print("KEY FEATURES DEMONSTRATED:")
print("=" * 60)
print("* 3D frame geometry definition")
print("* Nodal loads (wind force at node)")
print("* Member distributed loads (gravity on beam)")
print("* Member point loads (concentrated load at midspan)")
print("* Direct stiffness method analysis")
print("* Support reactions and equilibrium verification")
print("* Member force/stress/deflection results")
print("* 3D wireframe visualization with color mapping")

