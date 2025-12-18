"""
Custom Frame Structure
- Base frame: 250x150x6 RHS members forming perimeter
- 4 vertical posts: 150x150x6 RHS rising 1775mm
- 2 horizontal connecting posts: 150x150x6 RHS at 300mm height
"""

import numpy as np
from pathlib import Path
from sectiony.library import rhs, chs

# Import frame module
from beamy import Material
from beamy.frame import (
    Node, Member, Frame, FrameLoadCase, LoadedFrame
)

# Create gallery subdirectories
gallery_dir = Path("gallery")
gallery_dir.mkdir(exist_ok=True)
frame_dir = gallery_dir / "frame"
frame_dir.mkdir(exist_ok=True)

# ============================================
# 1. DEFINE MATERIALS AND SECTIONS
# ============================================
steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
sling = Material(name="Sling", E=50e9, G=20e9, Fy=1000e6)

# Base frame sections (250x150x6 RHS)
# Note: dimensions in mm, convert to m
base_section = rhs(b=0.150, h=0.250, t=0.006, r=0.0)

# Post sections (150x150x6 RHS)
post_section = rhs(b=0.150, h=0.150, t=0.006, r=0.0)

# Sling section (small CHS to approximate a cable)
sling_section = chs(d=0.020, t=0.002)

# ============================================
# 2. DEFINE FRAME GEOMETRY
# ============================================
# Convert mm to m for consistency
mm_to_m = 0.001

# Define all unique node positions
# Note: Members 13 & 14 connect to posts 9-12 at Z=300, so posts need nodes there
nodes = [
    # Member 1-3 nodes (left rail, X=125)
    Node("N1", np.array([125, 0, 0]) * mm_to_m),
    Node("N2", np.array([125, 750, 0]) * mm_to_m),
    Node("N3", np.array([125, 1750, 0]) * mm_to_m),
    Node("N4", np.array([125, 2500, 0]) * mm_to_m),
    
    # Member 4-6 nodes (right rail, X=1925)
    Node("N5", np.array([1925, 0, 0]) * mm_to_m),
    Node("N6", np.array([1925, 750, 0]) * mm_to_m),
    Node("N7", np.array([1925, 1750, 0]) * mm_to_m),
    Node("N8", np.array([1925, 2500, 0]) * mm_to_m),
    
    # Member 7 nodes (cross beam at Y=750)
    Node("N9", np.array([0, 750, 0]) * mm_to_m),
    Node("N10", np.array([2050, 750, 0]) * mm_to_m),
    
    # Member 8 nodes (cross beam at Y=1750)
    Node("N11", np.array([0, 1750, 0]) * mm_to_m),
    Node("N12", np.array([2050, 1750, 0]) * mm_to_m),
    
    # Post nodes (ground, 300mm, and top)
    # Post 1 (X=600, Y=750): Member 9 + part of Member 13
    Node("N13", np.array([600, 750, 0]) * mm_to_m),      # Ground
    Node("N21", np.array([600, 750, 300]) * mm_to_m),                      # 300mm (connection point)
    Node("N14", np.array([600, 750, 1775]) * mm_to_m),                     # Top
    
    # Post 2 (X=600, Y=1750): Member 10 + part of Member 13
    Node("N15", np.array([600, 1750, 0]) * mm_to_m),     # Ground
    Node("N22", np.array([600, 1750, 300]) * mm_to_m),                     # 300mm (connection point)
    Node("N16", np.array([600, 1750, 1775]) * mm_to_m),                    # Top
    
    # Post 3 (X=1450, Y=750): Member 11 + part of Member 14
    Node("N17", np.array([1450, 750, 0]) * mm_to_m),     # Ground
    Node("N23", np.array([1450, 750, 300]) * mm_to_m),                     # 300mm (connection point)
    Node("N18", np.array([1450, 750, 1775]) * mm_to_m),                    # Top
    
    # Post 4 (X=1450, Y=1750): Member 12 + part of Member 14
    Node("N19", np.array([1450, 1750, 0]) * mm_to_m),    # Ground
    Node("N24", np.array([1450, 1750, 300]) * mm_to_m),                    # 300mm (connection point)
    Node("N20", np.array([1450, 1750, 1775]) * mm_to_m),                   # Top

    # Lifting node (centered, 1200mm above column tops)
    Node("LIFT", np.array([1025, 1250, 1775 + 1200]) * mm_to_m),
]

# Define the 14 main structural members
# Note: Vertical posts (9-12) must be split into 2 segments each to allow connection at Z=300
members = [
    # Base frame (250x150x6 RHS)
    # Member.constraints is a 12-digit string like: [start 6][end 6]
    # where each digit corresponds to [UX UY UZ RX RY RZ] and 1=fixed.
    # Here we fully fix both ends of the base members to emulate ground anchorage.
    Member("M1", "N1", "N2", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (125,0,0) to (125,750,0)
    Member("M2", "N2", "N3", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (125,750,0) to (125,1750,0)
    Member("M3", "N3", "N4", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (125,1750,0) to (125,2500,0)
    Member("M4", "N5", "N6", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (1925,0,0) to (1925,750,0)
    Member("M5", "N6", "N7", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (1925,750,0) to (1925,1750,0)
    Member("M6", "N7", "N8", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),   # (1925,1750,0) to (1925,2500,0)
    Member("M7", "N9", "N10", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"),  # (0,750,0) to (2050,750,0)
    Member("M8", "N11", "N12", base_section, steel, np.array([0, 0, 1]), constraints="111111111111"), # (0,1750,0) to (2050,1750,0)
    
    # Vertical posts - LOWER SEGMENTS (ground to 300mm)
    # Fix only the ground end (start) of each post segment.
    Member("M9", "N13", "N21", post_section, steel, np.array([1, 0, 0]), constraints="111111000000"),  # Post 1 lower (600,750,0) to (600,750,300)
    Member("M10", "N15", "N22", post_section, steel, np.array([1, 0, 0]), constraints="111111000000"), # Post 2 lower (600,1750,0) to (600,1750,300)
    Member("M11", "N17", "N23", post_section, steel, np.array([1, 0, 0]), constraints="111111000000"), # Post 3 lower (1450,750,0) to (1450,750,300)
    Member("M12", "N19", "N24", post_section, steel, np.array([1, 0, 0]), constraints="111111000000"), # Post 4 lower (1450,1750,0) to (1450,1750,300)
    
    # Connecting posts at 300mm height (150x150x6 RHS)
    # Release all rotations (pinned ends) so crossbars don't restrain post bending
    # Format: 12 digits [start: UX UY UZ RX RY RZ][end: UX UY UZ RX RY RZ]
    # 0=rigid, 1=released. Release rotations: "000111000111"
    Member("M13", "N21", "N22", post_section, steel, np.array([0, 0, 1]), releases="000111000111"), # (600,750,300) to (600,1750,300)
    Member("M14", "N23", "N24", post_section, steel, np.array([0, 0, 1]), releases="000111000111"), # (1450,750,300) to (1450,1750,300)
    
    # Vertical posts - UPPER SEGMENTS (300mm to top)
    Member("M9_upper", "N21", "N14", post_section, steel, np.array([1, 0, 0])),  # Post 1 upper (600,750,300) to (600,750,1775)
    Member("M10_upper", "N22", "N16", post_section, steel, np.array([1, 0, 0])), # Post 2 upper (600,1750,300) to (600,1750,1775)
    Member("M11_upper", "N23", "N18", post_section, steel, np.array([1, 0, 0])), # Post 3 upper (1450,750,300) to (1450,750,1775)
    Member("M12_upper", "N24", "N20", post_section, steel, np.array([1, 0, 0])), # Post 4 upper (1450,1750,300) to (1450,1750,1775)

    # Sling members (4) from column tops to lifting node
    # releases string is 12 digits [start: UX UY UZ RX RY RZ][end: UX UY UZ RX RY RZ]
    # Release torsion only (RX) at both ends: "000100000100"
    Member("SL1", "N14", "LIFT", sling_section, sling, np.array([0, 0, 1]), releases="000100000100"),
    Member("SL2", "N16", "LIFT", sling_section, sling, np.array([0, 0, 1]), releases="000100000100"),
    Member("SL3", "N18", "LIFT", sling_section, sling, np.array([0, 0, 1]), releases="000100000100"),
    Member("SL4", "N20", "LIFT", sling_section, sling, np.array([0, 0, 1]), releases="000100000100"),
]

# Create frame
frame = Frame.from_nodes_and_members(nodes, members)

print("=" * 60)
print("CUSTOM FRAME STRUCTURE")
print("=" * 60)
print(f"Frame: {len(frame.nodes)} nodes, {len(frame.members)} members")
print(f"\nMember lengths:")
for member_id, length in frame.member_lengths.items():
    print(f"  {member_id}: {length*1000:.1f} mm ({length:.3f} m)")

# ============================================
# 3. DEFINE LOADS (Example: gravity on top posts)
# ============================================
loads = FrameLoadCase("Gravity + Wind")

# Example: Apply distributed gravity load on vertical posts
# Each post carries some weight (example: 1 kN/m)
gravity_per_meter = 1000.0  # 1 kN/m in N/m
print(f"\n{'=' * 60}")
print("APPLIED LOADS")
print("=" * 60)
print(f"Gravity load on vertical posts: {gravity_per_meter/1000:.1f} kN/m (global -Z)")

# Apply to all vertical post segments (lower and upper)
for post_id in ["M9", "M10", "M11", "M12", "M9_upper", "M10_upper", "M11_upper", "M12_upper"]:
    loads.add_member_uniform_force(post_id, np.array([0.0, 0.0, -gravity_per_meter]))

# Example: Add lateral wind load
wind_load = 5000.0  # 5 kN
print(f"Wind load at top nodes: {wind_load/1000:.1f} kN (global +X)")
for top_node in ["N14", "N16", "N18", "N20"]:  # Top of posts
    loads.add_nodal_force(top_node, np.array([wind_load, 0.0, 0.0]))

# ============================================
# 4. ANALYZE FRAME
# ============================================
print(f"\n{'=' * 60}")
print("PERFORMING FRAME ANALYSIS")
print("=" * 60)
print("Using direct stiffness method...")

loaded_frame = LoadedFrame(frame, loads)

print("Analysis complete!")

# ============================================
# 5. EXTRACT RESULTS
# ============================================
print(f"\n{'=' * 60}")
print("SUPPORT REACTIONS")
print("=" * 60)

total_reaction = np.zeros(3)
for node_id, reaction in loaded_frame.reactions.items():
    forces = reaction[:3]  # [FX, FY, FZ]
    moments = reaction[3:]  # [MX, MY, MZ]
    total_reaction += forces
    
    # Only print if significant reaction
    if np.linalg.norm(forces) > 1e-6:
        print(f"Node {node_id}:")
        print(f"  Forces:  FX={forces[0]:8.2f} N, FY={forces[1]:8.2f} N, FZ={forces[2]:8.2f} N")

print(f"\nTotal reaction forces:")
print(f"  FX={total_reaction[0]:8.2f} N")
print(f"  FY={total_reaction[1]:8.2f} N")
print(f"  FZ={total_reaction[2]:8.2f} N")

# ============================================
# 6. VISUALIZE
# ============================================
print(f"\n{'=' * 60}")
print("GENERATING VISUALIZATION PLOTS")
print("=" * 60)

# Basic geometry with loads
print("1. Frame geometry...")
loaded_frame.plot(
    deformed=True,
    scale_factor=100,
    save_path=str(frame_dir / "custom_frame_geometry.svg")
)

# Deflection plot
print("2. Deflection plot...")
loaded_frame.plot_deflection(
    scale_factor=100,
    colormap="viridis",
    show_undeformed=True,
    save_path=str(frame_dir / "custom_frame_deflection.svg")
)

# Von Mises stress plot
print("3. Von Mises stress distribution...")
loaded_frame.plot_von_mises(
    colormap="turbo",
    stress_limits=(0, steel.Fy),
    save_path=str(frame_dir / "custom_frame_von_mises.svg")
)

print(f"\n{'=' * 60}")
print("COMPLETE!")
print("=" * 60)
print(f"\nPlots saved to: {frame_dir.absolute()}/")
print("  - custom_frame_geometry.svg")
print("  - custom_frame_deflection.svg")
print("  - custom_frame_von_mises.svg")


