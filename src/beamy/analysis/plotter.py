from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

if TYPE_CHECKING:
    from ..setup import Beam1D, LoadCase


def _plot_moments(ax, moments):
    pass

def _plot_distributed_forces(ax, dist_forces):
    pass

def _create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=8):
    """
    Creates a 3D cone for the arrow head.
    
    Args:
        tip: 3D point where the cone tip is located
        direction: Normalized direction vector (points from base to tip)
        cone_length: Length of the cone
        cone_radius: Radius of the cone base
        num_segments: Number of segments to approximate the cone base
    
    Returns:
        List of triangular faces for the cone
    """
    # Base center is at tip - direction * cone_length
    base_center = tip - direction * cone_length
    
    # Create two perpendicular vectors to the direction for the base circle
    # Find a vector perpendicular to direction
    if abs(direction[0]) < 0.9:
        perp1 = np.array([1, 0, 0])
    else:
        perp1 = np.array([0, 1, 0])
    perp1 = perp1 - np.dot(perp1, direction) * direction
    perp1 = perp1 / np.linalg.norm(perp1)
    
    # Second perpendicular vector
    perp2 = np.cross(direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # Generate base circle points
    base_points = []
    for i in range(num_segments):
        angle = 2 * np.pi * i / num_segments
        point = base_center + cone_radius * (np.cos(angle) * perp1 + np.sin(angle) * perp2)
        base_points.append(point)
    
    # Create triangular faces (each triangle: tip, base[i], base[i+1])
    faces = []
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        faces.append([tip, base_points[i], base_points[next_i]])
    
    return faces

def _plot_point_forces(ax, point_forces, beam_length):
    """
    Plots 3D arrows for point forces using lines and 3D cones.
    Mapping:
    Beam X -> Plot Y
    Beam Y -> Plot Z
    Beam Z -> Plot X
    
    Args:
        ax: 3D axes object
        point_forces: List of PointForce objects
        beam_length: Length of the beam (used for scaling arrow lengths)
    """
    if not point_forces:
        return
    
    # Calculate force magnitudes to determine scaling
    magnitudes = [np.linalg.norm(pf.force) for pf in point_forces]
    max_magnitude = max(magnitudes) if magnitudes else 1.0
    
    # Longest arrow should be 1/5th the length of the beam
    max_arrow_length = beam_length / 5.0
    
    for pf in point_forces:
        # Extract Beam coordinates and forces
        # PointForce.point is [x, y, z]
        px, py, pz = pf.point
        # PointForce.force is [Fx, Fy, Fz]
        fx, fy, fz = pf.force
        
        # Map to Plot coordinates
        # Plot (X, Y, Z) = (Beam Z, Beam X, Beam Y)
        tip = np.array([pz, px, py])  # Tip is at the point
        force_vec = np.array([fz, fx, fy])  # Force direction in plot coordinates
        
        # Calculate arrow length proportional to force magnitude
        force_magnitude = np.linalg.norm(force_vec)
        if max_magnitude > 0:
            arrow_length = (force_magnitude / max_magnitude) * max_arrow_length
        else:
            arrow_length = max_arrow_length
        
        # Normalize direction vector
        if force_magnitude > 0:
            direction = force_vec / force_magnitude
        else:
            direction = np.array([1, 0, 0])  # Default direction
        
        # Tail is at tip - direction * arrow_length (since tip is at point)
        tail = tip - direction * arrow_length
        
        # Draw the arrow shaft (line)
        ax.plot([tail[0], tip[0]], 
                [tail[1], tip[1]], 
                [tail[2], tip[2]], 
                color='red', linewidth=2)
        
        # Create and draw the cone (arrow head)
        cone_length = arrow_length * 0.15  # Cone is 15% of arrow length
        cone_radius = arrow_length * 0.05   # Cone radius is 5% of arrow length
        cone_faces = _create_arrow_cone(tip, direction, cone_length, cone_radius)
        
        cone_collection = Poly3DCollection(cone_faces, facecolor='red', edgecolor='red', alpha=1.0)
        ax.add_collection3d(cone_collection)

def plot_loads(ax, load_case: "LoadCase", beam_length: float):
    if not load_case:
        return
    
    _plot_point_forces(ax, load_case.point_forces, beam_length)
    _plot_moments(ax, load_case.moments)
    _plot_distributed_forces(ax, load_case.dist_forces)

def plot_beam_diagram(beam: "Beam1D", load_case: Optional["LoadCase"] = None):
    """
    Plots a 3D extrusion of the beam section.
    The section is defined in the YZ plane and extruded along the X axis.
    
    Args:
        beam: Beam1D object containing section geometry and length.
        load_case: Optional LoadCase object containing forces to plot.
    """
    if not beam.section.geometry or not beam.section.geometry.shapes:
        print("Warning: Beam section has no geometry defined. Cannot plot 3D extrusion.")
        return

    # For now, we only plot the first shape in the geometry (ignoring hollows/multi-shape for basic plot)
    # Ideally, we should iterate over all shapes or handle complex geometry.
    # Assuming first shape is the outer boundary for simple visualization.
    shape_points = beam.section.geometry.shapes[0].points
    length = beam.L
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Convert to numpy array for easier handling
    pts = np.array(shape_points)
    
    # Validate inputs
    if len(pts) < 3:
        print("Need at least 3 points to define a shape")
        return

    # Vertices at x=0 (Start of beam)
    # Mapping to plot coordinates (User request):
    # Plot X <= Beam Z
    # Plot Y <= Beam X
    # Plot Z <= Beam Y
    x0 = np.zeros(pts.shape[0])
    # (BeamZ, BeamX, BeamY)
    base_verts = np.column_stack((pts[:, 1], x0, pts[:, 0]))

    # Vertices at x=length (End of beam)
    x_end = np.full(pts.shape[0], length)
    end_verts = np.column_stack((pts[:, 1], x_end, pts[:, 0]))

    faces = []

    # 1. Start face (Beam x=0 -> Plot Y=0)
    faces.append(base_verts)

    # 2. End face (Beam x=length -> Plot Y=length)
    faces.append(end_verts)

    # 3. Side faces
    # Connect corresponding points between base and end
    # Iterate through segments. Since points are closed, we have len-1 segments.
    for i in range(len(pts) - 1):
        # Define the quad for the side face
        # Order: Base[i], Base[i+1], End[i+1], End[i]
        side = [
            base_verts[i],
            base_verts[i+1],
            end_verts[i+1],
            end_verts[i]
        ]
        faces.append(side)
    
    # If the loop is closed (first point == last point), loop logic covers it.
    # If not explicitly closed in data, we might miss the last closing face.
    # Assuming standard polygon definition where first != last usually implies implicit closure, 
    # but here we iterate `range(len(pts) - 1)` which implies pts has the start point repeated at end.
    # If not repeated, we need to close it manually.
    if not np.allclose(pts[0], pts[-1]):
        # Close the loop
        side = [
            base_verts[-1],
            base_verts[0],
            end_verts[0],
            end_verts[-1]
        ]
        faces.append(side)

    # Create the collection
    mesh = Poly3DCollection(faces, alpha=0.5, edgecolor='k')
    
    # Optional: Set face colors
    mesh.set_facecolor('silver')

    ax.add_collection3d(mesh)
    
    # Plot Loads if provided
    if load_case:
        plot_loads(ax, load_case, length)

    # Set plot limits
    # We need to manually set limits because add_collection3d doesn't auto-scale
    all_coords = np.vstack((base_verts, end_verts))
    
    x_min, x_max = all_coords[:, 0].min(), all_coords[:, 0].max()
    y_min, y_max = all_coords[:, 1].min(), all_coords[:, 1].max()
    z_min, z_max = all_coords[:, 2].min(), all_coords[:, 2].max()
    
    # Add some padding
    max_range = np.array([x_max-x_min, y_max-y_min, z_max-z_min]).max() / 2.0
    mid_x = (x_max+x_min) * 0.5
    mid_y = (y_max+y_min) * 0.5
    mid_z = (z_max+z_min) * 0.5
    
    # Ensure limits are at least some minimal size to avoid singular matrix errors
    if max_range == 0: max_range = 1.0

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('Z (Transverse Horizontal)')
    ax.set_ylabel('X (Longitudinal)')
    ax.set_zlabel('Y (Vertical)')
    
    # Set white background
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('white')
    ax.yaxis.pane.set_edgecolor('white')
    ax.zaxis.pane.set_edgecolor('white')
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5)
    
    ax.view_init(elev=20, azim=-60)
    
    plt.title(f"Extruded Section (L={length})")
    plt.show()
