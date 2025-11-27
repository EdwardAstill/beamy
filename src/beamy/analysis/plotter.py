from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import to_rgba, Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np

if TYPE_CHECKING:
    from .analysis import LoadedBeam


def _create_moment_arc(center, moment_vec_beam, radius, arc_angle=3*np.pi/2, num_points=30):
    """
    Creates a 3/4 arc in a plane perpendicular to the moment vector.
    
    Args:
        center: 3D point at the center of the arc (in plot coordinates)
        moment_vec_beam: Moment vector in BEAM coordinates [Mx, My, Mz]
        radius: Radius of the arc
        arc_angle: Total arc angle (default 3/4 turn = 270 degrees)
        num_points: Number of points to approximate the arc
    
    Returns:
        Tuple of (arc_points, end_direction) where end_direction is tangent at arc end
    """
    # Map moment vector to plot coordinates
    # Beam (Mx, My, Mz) -> Plot (Mz, Mx, My)
    m_plot = np.array([moment_vec_beam[2], moment_vec_beam[0], moment_vec_beam[1]])
    
    # Normalize normal vector
    norm = np.linalg.norm(m_plot)
    if norm < 1e-10:
        return np.array([]), np.array([0,0,0])
    
    n = m_plot / norm
    
    # Find two orthogonal vectors u, v in the plane perpendicular to n
    # If n is roughly parallel to [1,0,0], use [0,1,0] as temp
    if abs(n[0]) < 0.9:
        u = np.array([1, 0, 0])
    else:
        u = np.array([0, 1, 0])
        
    u = u - np.dot(u, n) * n
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    
    # Create arc points
    # We want to spiral "around" the vector n according to right-hand rule
    # Rotation from u to v is rotation around n
    angles = np.linspace(0, arc_angle, num_points)
    
    arc_points = []
    for a in angles:
        # Point on circle in plane: center + R * (cos(a)*u + sin(a)*v)
        pt = center + radius * (np.cos(a) * u + np.sin(a) * v)
        arc_points.append(pt)
    
    arc_points = np.array(arc_points)
    
    # Tangent direction at end: -sin(end)*u + cos(end)*v
    # This is the direction of motion along the arc
    end_a = arc_angle
    end_dir = -np.sin(end_a) * u + np.cos(end_a) * v
    
    return arc_points, end_dir


def _plot_moments(ax, moments, beam_length):
    """
    Plots 3D moment arrows as 3/4 arcs with cone tips.
    
    Args:
        ax: 3D axes object
        moments: List of Moment objects
        beam_length: Length of the beam (used for scaling)
    """
    if not moments:
        return
    
    # Calculate max moment magnitude for scaling
    magnitudes = [np.linalg.norm(m.moment) for m in moments]
    max_magnitude = max(magnitudes) if magnitudes else 1.0
    
    if max_magnitude == 0:
        return
    
    # Max arc radius is 1/10 of beam length
    max_arc_radius = beam_length / 10.0
    
    for m in moments:
        magnitude = np.linalg.norm(m.moment)
        if magnitude < 1e-10:
            continue
            
        # Position along beam: m.x -> Plot Y coordinate
        # Plot coords: (Beam Z, Beam X, Beam Y) = (0, m.x, 0) for on-axis
        center = np.array([0, m.x, 0])
        
        # Scale radius by magnitude
        radius = (magnitude / max_magnitude) * max_arc_radius
        
        # Create arc (one arc per moment vector)
        # Pass the full moment vector
        arc_points, end_dir = _create_moment_arc(
            center, m.moment, radius, arc_angle=3*np.pi/2
        )
        
        if len(arc_points) == 0:
            continue
        
        # Draw cone at end of arc
        cone_length = radius * 0.3
        cone_radius = radius * 0.12
        
        # Tip is at arc end
        arc_end = arc_points[-1]
        
        # Shift cone so base center is at arc end
        tip = arc_end + end_dir * cone_length
        
        side_faces, base_face = _create_arrow_cone(
            tip, end_dir, cone_length, cone_radius, num_segments=12
        )
        
        # Draw in order: side faces, then base, then arc
        # Side faces in blue
        side_collection = Poly3DCollection(side_faces, facecolor='blue', edgecolor='blue', alpha=1.0)
        ax.add_collection3d(side_collection)
        
        # Base face in lighter blue
        base_collection = Poly3DCollection(base_face, facecolor='#4444dd', edgecolor='#4444dd', alpha=1.0)
        ax.add_collection3d(base_collection)
        
        # Draw arc (on top)
        ax.plot(arc_points[:, 0], arc_points[:, 1], arc_points[:, 2],
                color='blue', linewidth=2, zorder=100)


def _plot_distributed_forces(ax, dist_forces, beam_length):
    """
    Plots 3D arrows for distributed forces using lines and 3D cones.
    Draws arrows at start and end, plus intermediate arrows based on load length ratio.
    Connects all arrow tails with a line.
    
    Mapping:
    Beam X -> Plot Y
    Beam Y -> Plot Z
    Beam Z -> Plot X
    
    Args:
        ax: 3D axes object
        dist_forces: List of DistributedForce objects
        beam_length: Length of the beam (used for scaling arrow lengths)
    """
    if not dist_forces:
        return
    
    # Calculate max force magnitude across all distributed forces for scaling
    max_magnitude = 0.0
    for df in dist_forces:
        max_magnitude = max(max_magnitude, np.linalg.norm(df.start_force), np.linalg.norm(df.end_force))
    
    if max_magnitude == 0:
        return
    
    # Longest arrow should be 1/5th the length of the beam
    max_arrow_length = beam_length / 5.0
    
    for df in dist_forces:
        # Calculate distributed load length
        load_length = np.linalg.norm(df.end_position - df.start_position)
        if load_length < 1e-10:
            continue
        
        # Calculate average force magnitude from start and end for consistent cone size
        start_force_mag = np.linalg.norm(df.start_force)
        end_force_mag = np.linalg.norm(df.end_force)
        avg_force_mag = (start_force_mag + end_force_mag) / 2.0
        
        if avg_force_mag < 1e-10:
            continue
        
        # Calculate average arrow length based on average force magnitude
        avg_arrow_length = (avg_force_mag / max_magnitude) * max_arrow_length
        
        # Fixed cone size for all arrows in this distributed load
        cone_length = avg_arrow_length * 0.2
        cone_radius = avg_arrow_length * 0.08
        
        # Determine number of arrows based on ratio of load length to beam length
        # Base count of 10 arrows for full beam length, minimum 2 (start and end)
        base_count = 10
        num_arrows = max(2, int(np.ceil((load_length / beam_length) * base_count)))
        
        # Interpolate positions and forces
        ts = np.linspace(0.0, 1.0, num_arrows)
        positions = np.outer(1 - ts, df.start_position) + np.outer(ts, df.end_position)
        forces = np.outer(1 - ts, df.start_force) + np.outer(ts, df.end_force)
        
        # Store tail positions for connecting line
        tail_positions = []
        
        for i in range(num_arrows):
            pos = positions[i]
            force = forces[i]
            
            # Extract Beam coordinates
            px, py, pz = pos
            fx, fy, fz = force
            
            # Map to Plot coordinates
            # Plot (X, Y, Z) = (Beam Z, Beam X, Beam Y)
            tip = np.array([pz, px, py])
            force_vec = np.array([fz, fx, fy])
            
            # Calculate arrow length proportional to force magnitude
            force_magnitude = np.linalg.norm(force_vec)
            if force_magnitude < 1e-10:
                continue
            
            arrow_length = (force_magnitude / max_magnitude) * max_arrow_length
            
            # Normalize direction vector
            direction = force_vec / force_magnitude
            
            # Tail is at tip - direction * arrow_length
            tail = tip - direction * arrow_length
            tail_positions.append(tail)
            
            # Create and draw the cone (arrow head) with fixed size
            side_faces, base_face = _create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=12)
            
            # Side faces in green
            side_collection = Poly3DCollection(side_faces, facecolor='green', edgecolor='green', alpha=1.0)
            ax.add_collection3d(side_collection)
            
            # Base face in lighter green
            base_collection = Poly3DCollection(base_face, facecolor='#44dd44', edgecolor='#44dd44', alpha=1.0)
            ax.add_collection3d(base_collection)
            
            # Draw the arrow shaft
            cone_base = tip - direction * cone_length
            ax.plot([tail[0], cone_base[0]], 
                    [tail[1], cone_base[1]], 
                    [tail[2], cone_base[2]], 
                    color='green', linewidth=2, zorder=100)
        
        # Connect all arrow tails with a line
        if len(tail_positions) >= 2:
            tail_positions = np.array(tail_positions)
            ax.plot(tail_positions[:, 0], tail_positions[:, 1], tail_positions[:, 2],
                    color='green', linewidth=2, zorder=100)


def _plot_supports(ax, supports, beam_length):
    """
    Plots supports as small spheres with their 6-digit code labels.
    
    Args:
        ax: 3D axes object
        supports: List of Support objects
        beam_length: Length of the beam (used for scaling sphere size)
    """
    if not supports:
        return
    
    # Sphere radius scales with beam length
    sphere_radius = beam_length / 100.0
    
    for support in supports:
        # Support position along beam axis
        # Plot coords: (Beam Z, Beam X, Beam Y) = (0, support.x, 0)
        center = np.array([0, support.x, 0])
        
        # Create sphere mesh
        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 10)
        x = center[0] + sphere_radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + sphere_radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + sphere_radius * np.outer(np.ones(np.size(u)), np.cos(v))
        
        # Plot sphere surface
        ax.plot_surface(x, y, z, color='black', alpha=0.8, linewidth=0)
        
        # Add label with 6-digit code, offset slightly above the sphere
        label_offset = sphere_radius * 1.5
        ax.text(center[0], center[1], center[2] + label_offset, 
                support.type, fontsize=8, ha='center', va='bottom', color='black')

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
        Tuple of (side_faces, base_face) - side triangles and base polygon
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
    
    # Create triangular side faces (each triangle: tip, base[i], base[i+1])
    side_faces = []
    for i in range(num_segments):
        next_i = (i + 1) % num_segments
        side_faces.append([tip, base_points[i], base_points[next_i]])
    
    # Base face (polygon of all base points)
    base_face = [base_points]
    
    return side_faces, base_face

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
        
        # Create and draw the cone (arrow head) FIRST
        cone_length = arrow_length * 0.2   # Cone is 20% of arrow length
        cone_radius = arrow_length * 0.08  # Cone radius is 8% of arrow length
        side_faces, base_face = _create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=12)
        
        cone_base = tip - direction * cone_length
        
        # Side faces in red
        side_collection = Poly3DCollection(side_faces, facecolor='red', edgecolor='red', alpha=1.0)
        ax.add_collection3d(side_collection)
        
        # Base face
        base_collection = Poly3DCollection(base_face, facecolor='#ff4444', edgecolor='#ff4444', alpha=1.0)
        ax.add_collection3d(base_collection)
        
        # Draw the arrow shaft line last so it appears on top
        ax.plot([tail[0], cone_base[0]], 
                [tail[1], cone_base[1]], 
                [tail[2], cone_base[2]], 
                color='red', linewidth=2, zorder=100)

def plot_loads(ax, load_case, beam_length: float):
    if not load_case:
        return
    
    _plot_point_forces(ax, load_case.point_forces, beam_length)
    _plot_moments(ax, load_case.moments, beam_length)
    _plot_distributed_forces(ax, load_case.dist_forces, beam_length)


def _plot_stress_line(ax, loaded_beam: "LoadedBeam", length: float, n_points: int = 100):
    """
    Plot the beam axis as a color-coded line based on von Mises stress.
    
    Uses conservative max von Mises stress along the beam (including shear):
    σ_vm = √(σ² + 3τ²)
    """
    # Get von Mises stress results directly from analysis
    # This includes axial, bending, shear, and torsion
    vm_result = loaded_beam.von_mises(points=n_points)
    
    xs = vm_result._x
    von_mises = vm_result._values
    
    # Create line segments for color mapping
    # Plot coords: (0, x, 0) for each beam x position
    points = np.array([[0, x, 0] for x in xs])
    segments = np.array([[points[i], points[i+1]] for i in range(len(points)-1)])
    
    # Color mapping
    norm = Normalize(vmin=von_mises.min(), vmax=von_mises.max())
    cmap = cm.plasma  # Blue (low) to red/yellow (high)
    
    # Use midpoint stress for each segment
    segment_stress = (von_mises[:-1] + von_mises[1:]) / 2
    colors = cmap(norm(segment_stress))
    
    # Create line collection
    lc = Line3DCollection(segments, colors=colors, linewidth=3)
    ax.add_collection3d(lc)
    
    # Add colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='Max Von Mises Stress', shrink=0.6, pad=0.1)

def plot_beam_diagram(loaded_beam: "LoadedBeam", plot_stress: bool = False, plot_section: bool = True):
    """
    Plots a 3D beam diagram with:
    - Optional 2D section outline on the YZ plane at x=0 (shear center at origin)
    - A line along the beam x-axis representing the beam length
    - Optional von Mises stress coloring on the beam axis
    
    Args:
        loaded_beam: LoadedBeam object containing beam, loads, and analysis results.
        plot_stress: If True, color the beam axis by von Mises stress.
        plot_section: If True, draw the section outline at x=0.
    """
    beam = loaded_beam.beam
    load_case = loaded_beam.loads
    
    if not beam.section.geometry or not beam.section.geometry.shapes:
        print("Warning: Beam section has no geometry defined. Cannot plot.")
        return

    length = beam.L
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Get shear center offsets (to position section with shear center at origin)
    sc_y = getattr(beam.section, 'SCy', 0.0)
    sc_z = getattr(beam.section, 'SCz', 0.0)

    # Draw 2D section outline on the YZ plane at beam x=0
    if plot_section:
        # Coordinate mapping: Plot (X, Y, Z) = (Beam Z, Beam X, Beam Y)
        for shape in beam.section.geometry.shapes:
            pts = np.array(shape.points)
            if len(pts) < 3:
                continue
            
            # Offset points so shear center is at origin
            # pts[:, 0] is beam Y, pts[:, 1] is beam Z
            offset_y = pts[:, 0] - sc_y
            offset_z = pts[:, 1] - sc_z
            
            # Close the polygon if not already closed
            if not np.allclose(pts[0], pts[-1]):
                offset_y = np.append(offset_y, offset_y[0])
                offset_z = np.append(offset_z, offset_z[0])
            
            # Plot section at beam x=0 -> Plot Y=0
            # Plot coords: (Beam Z, Beam X, Beam Y) = (offset_z, 0, offset_y)
            plot_x = offset_z
            plot_y = np.zeros_like(offset_z)
            plot_z = offset_y
            
            ax.plot(plot_x, plot_y, plot_z, color='grey', linewidth=1.5)

    # Draw beam axis line from (0,0,0) to (L,0,0) in beam coords
    # Plot coords: (0, 0, 0) to (0, L, 0)
    if plot_stress:
        _plot_stress_line(ax, loaded_beam, length)
    else:
        ax.plot([0, 0], [0, length], [0, 0], color='black', linewidth=2)
    
    # Plot Loads
    plot_loads(ax, load_case, length)

    # Set plot limits based on section size and beam length
    # Get max section extent for scaling
    all_y = []
    all_z = []
    for shape in beam.section.geometry.shapes:
        pts = np.array(shape.points)
        all_y.extend(pts[:, 0] - sc_y)
        all_z.extend(pts[:, 1] - sc_z)
    
    section_extent = max(abs(min(all_y)), abs(max(all_y)), abs(min(all_z)), abs(max(all_z)))
    if section_extent == 0:
        section_extent = 1.0
    
    # Use the larger of section extent or fraction of beam length for YZ scaling
    yz_range = max(section_extent * 1.5, length * 0.1)
    
    ax.set_xlim(-yz_range, yz_range)
    ax.set_ylim(0, length)
    ax.set_zlim(-yz_range, yz_range)
    
    # Set equal aspect ratio so cones appear as proper cones
    ax.set_box_aspect([2 * yz_range, length, 2 * yz_range])

    # Hide all panes and grid
    ax.grid(False)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    # Hide Plot X axis completely (Beam Z)
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.set_xticks([])
    
    # Plot Y axis = Beam X (along beam) - label as "x", show ticks, position further from axis
    ax.set_ylabel('x', labelpad=20)
    
    # Plot Z axis = Beam Y (vertical) - label as "y", hide ticks, position closer to axis
    ax.set_zlabel('y', labelpad=-10)
    ax.set_zticks([])
    
    # Plot Supports last so they overlap everything
    _plot_supports(ax, beam.supports, length)
    
    # Set white background
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    ax.view_init(elev=20, azim=-60)
    
    title = f"Beam Diagram (L={length})"
    if plot_stress:
        title = f"Von Mises Stress (L={length})"
    plt.title(title)
    plt.show()
