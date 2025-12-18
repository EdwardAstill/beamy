# plotter.py
from __future__ import annotations
from typing import Optional, Tuple, TYPE_CHECKING
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

if TYPE_CHECKING:
    from .analysis import LoadedFrame


def plot_frame(
    loaded_frame: 'LoadedFrame',
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
        loaded_frame: LoadedFrame object with analysis results
        show_loads: Display applied load arrows
        show_reactions: Display reaction arrows at supports
        show_member_ids: Label members with their IDs
        show_node_ids: Label nodes with their IDs
        deformed: If True, show deformed shape overlay
        scale_factor: Scale factor for deformed shape visualization
        save_path: Path to save the plot (supports .svg)
    """
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get node positions
    node_positions = loaded_frame.frame.node_positions
    
    # Plot undeformed frame
    for member in loaded_frame.frame.members:
        start_pos = node_positions[member.start_node_id]
        end_pos = node_positions[member.end_node_id]
        
        xs = [start_pos[0], end_pos[0]]
        ys = [start_pos[1], end_pos[1]]
        zs = [start_pos[2], end_pos[2]]
        
        ax.plot(xs, ys, zs, 'k-', linewidth=2, label='Undeformed' if member == loaded_frame.frame.members[0] else '')
    
    # Plot deformed shape if requested
    if deformed:
        for member in loaded_frame.frame.members:
            start_pos = node_positions[member.start_node_id]
            end_pos = node_positions[member.end_node_id]
            
            start_disp = loaded_frame.nodal_displacements[member.start_node_id][:3]
            end_disp = loaded_frame.nodal_displacements[member.end_node_id][:3]
            
            start_def = start_pos + scale_factor * start_disp
            end_def = end_pos + scale_factor * end_disp
            
            xs = [start_def[0], end_def[0]]
            ys = [start_def[1], end_def[1]]
            zs = [start_def[2], end_def[2]]
            
            ax.plot(xs, ys, zs, 'b--', linewidth=1.5, alpha=0.7, 
                   label='Deformed' if member == loaded_frame.frame.members[0] else '')
    
    # Plot nodes
    for node_id, pos in node_positions.items():
        node = loaded_frame.frame.nodes[node_id]
        if node.support is not None:
            # Supported nodes as triangles
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='red', marker='^', s=100, zorder=5)
        else:
            # Free nodes as circles
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='black', marker='o', s=50, zorder=5)
        
        # Node labels
        if show_node_ids:
            ax.text(pos[0], pos[1], pos[2], f'  {node_id}', fontsize=10, weight='bold')
    
    # Member labels
    if show_member_ids:
        for member in loaded_frame.frame.members:
            start_pos = node_positions[member.start_node_id]
            end_pos = node_positions[member.end_node_id]
            mid_pos = (start_pos + end_pos) / 2
            ax.text(mid_pos[0], mid_pos[1], mid_pos[2], f' {member.id}', 
                   fontsize=9, color='blue', style='italic')
    
    # Plot loads
    if show_loads:
        arrow_scale = _compute_arrow_scale(loaded_frame)
        
        # Nodal forces
        for nf in loaded_frame.loads.nodal_forces:
            pos = node_positions[nf.node_id]
            force_mag = np.linalg.norm(nf.force)
            if force_mag > 1e-10:
                force_dir = nf.force / force_mag
                ax.quiver(pos[0], pos[1], pos[2],
                         force_dir[0], force_dir[1], force_dir[2],
                         length=arrow_scale, color='red', arrow_length_ratio=0.3, linewidth=2)
    
    # Plot reactions
    if show_reactions:
        arrow_scale = _compute_arrow_scale(loaded_frame)
        
        for node_id, reaction in loaded_frame.reactions.items():
            pos = node_positions[node_id]
            force = reaction[:3]  # [FX, FY, FZ]
            force_mag = np.linalg.norm(force)
            if force_mag > 1e-10:
                force_dir = force / force_mag
                ax.quiver(pos[0], pos[1], pos[2],
                         force_dir[0], force_dir[1], force_dir[2],
                         length=arrow_scale, color='green', arrow_length_ratio=0.3, linewidth=2)
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Frame Structure')
    
    # Equal aspect ratio
    _set_axes_equal(ax)
    
    # Legend
    if deformed:
        ax.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    else:
        plt.show()


def plot_deflection(
    loaded_frame: 'LoadedFrame',
    scale_factor: float = 1.0,
    points_per_member: int = 20,
    colormap: str = "viridis",
    show_undeformed: bool = True,
    show_colorbar: bool = True,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the deformed frame shape in 3D wireframe, colored by displacement magnitude.
    
    Args:
        loaded_frame: LoadedFrame object with analysis results
        scale_factor: Multiplier for displacement visualization
        points_per_member: Number of interpolation points along each member
        colormap: Matplotlib colormap name
        show_undeformed: Show original geometry as faint dashed lines
        show_colorbar: Display colorbar with displacement units
        save_path: Path to save the plot (supports .svg)
    """
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    node_positions = loaded_frame.frame.node_positions
    
    # Collect all displacement magnitudes for color mapping
    all_disps = []
    
    for member in loaded_frame.frame.members:
        start_pos = node_positions[member.start_node_id]
        end_pos = node_positions[member.end_node_id]
        
        start_disp = loaded_frame.nodal_displacements[member.start_node_id][:3]
        end_disp = loaded_frame.nodal_displacements[member.end_node_id][:3]
        
        # Interpolate positions and displacements along member
        t = np.linspace(0, 1, points_per_member)
        for ti in t:
            interp_disp = (1 - ti) * start_disp + ti * end_disp
            disp_mag = np.linalg.norm(interp_disp)
            all_disps.append(disp_mag)
    
    # Setup color normalization
    vmin, vmax = min(all_disps), max(all_disps)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(colormap)
    
    # Plot undeformed if requested
    if show_undeformed:
        for member in loaded_frame.frame.members:
            start_pos = node_positions[member.start_node_id]
            end_pos = node_positions[member.end_node_id]
            
            xs = [start_pos[0], end_pos[0]]
            ys = [start_pos[1], end_pos[1]]
            zs = [start_pos[2], end_pos[2]]
            
            ax.plot(xs, ys, zs, 'k--', linewidth=1, alpha=0.3)
    
    # Plot deformed shape with color mapping
    for member in loaded_frame.frame.members:
        start_pos = node_positions[member.start_node_id]
        end_pos = node_positions[member.end_node_id]
        
        start_disp = loaded_frame.nodal_displacements[member.start_node_id][:3]
        end_disp = loaded_frame.nodal_displacements[member.end_node_id][:3]
        
        # Create line segments
        t = np.linspace(0, 1, points_per_member)
        points = []
        colors = []
        
        for ti in t:
            # Interpolate position and displacement
            interp_pos = (1 - ti) * start_pos + ti * end_pos
            interp_disp = (1 - ti) * start_disp + ti * end_disp
            
            # Deformed position
            def_pos = interp_pos + scale_factor * interp_disp
            points.append(def_pos)
            
            # Color by displacement magnitude
            disp_mag = np.linalg.norm(interp_disp)
            colors.append(cmap(norm(disp_mag)))
        
        # Plot as connected line segments with colors
        for i in range(len(points) - 1):
            xs = [points[i][0], points[i+1][0]]
            ys = [points[i][1], points[i+1][1]]
            zs = [points[i][2], points[i+1][2]]
            ax.plot(xs, ys, zs, color=colors[i], linewidth=2)
    
    # Add colorbar
    if show_colorbar:
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=10)
        cbar.set_label('Displacement Magnitude', rotation=270, labelpad=20)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Deflection (scale factor = {scale_factor})')
    
    _set_axes_equal(ax)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    else:
        plt.show()


def plot_von_mises(
    loaded_frame: 'LoadedFrame',
    points_per_member: int = 20,
    colormap: str = "jet",
    show_colorbar: bool = True,
    stress_limits: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> None:
    """
    Plot the frame in 3D wireframe, colored by Von Mises stress.
    
    Args:
        loaded_frame: LoadedFrame object with analysis results
        points_per_member: Number of interpolation points along each member
        colormap: Matplotlib colormap name
        show_colorbar: Display colorbar with stress units
        stress_limits: Optional (min, max) to fix colorbar range
        save_path: Path to save the plot (supports .svg)
    """
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    node_positions = loaded_frame.frame.node_positions
    
    # Collect all von Mises stresses
    all_stresses = []
    member_stress_data = {}
    
    for member in loaded_frame.frame.members:
        member_results = loaded_frame.get_member_results(member.id)
        vm_result = member_results.von_mises
        
        # Sample von Mises stress at points
        x_positions = np.linspace(0, member.length, points_per_member)
        stresses = [vm_result.at(x) for x in x_positions]
        
        member_stress_data[member.id] = stresses
        all_stresses.extend(stresses)
    
    # Setup color normalization
    if stress_limits:
        vmin, vmax = stress_limits
    else:
        vmin, vmax = min(all_stresses), max(all_stresses)
    
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(colormap)
    
    # Plot members with stress coloring
    for member in loaded_frame.frame.members:
        start_pos = node_positions[member.start_node_id]
        end_pos = node_positions[member.end_node_id]
        
        stresses = member_stress_data[member.id]
        t = np.linspace(0, 1, points_per_member)
        
        # Create line segments
        for i in range(len(t) - 1):
            t1, t2 = t[i], t[i+1]
            pos1 = (1 - t1) * start_pos + t1 * end_pos
            pos2 = (1 - t2) * start_pos + t2 * end_pos
            
            # Average stress for this segment
            stress_avg = (stresses[i] + stresses[i+1]) / 2
            color = cmap(norm(stress_avg))
            
            xs = [pos1[0], pos2[0]]
            ys = [pos1[1], pos2[1]]
            zs = [pos1[2], pos2[2]]
            
            ax.plot(xs, ys, zs, color=color, linewidth=3)
    
    # Add colorbar
    if show_colorbar:
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=10)
        cbar.set_label('Von Mises Stress (Pa)', rotation=270, labelpad=20)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Von Mises Stress Distribution')
    
    _set_axes_equal(ax)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    else:
        plt.show()


def plot_results(
    loaded_frame: 'LoadedFrame',
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
        loaded_frame: LoadedFrame object with analysis results
        result_type: Type of result to display
        deformed: Plot on deformed geometry
        scale_factor: Displacement scale factor
        points_per_member: Number of interpolation points per member
        colormap: Matplotlib colormap name
        show_undeformed: Show original geometry as reference
        show_colorbar: Display colorbar legend
        show_node_ids: Label nodes
        show_member_ids: Label members
        value_limits: Optional (min, max) to fix color range
        save_path: Path to save the plot
    """
    if result_type == "von_mises":
        plot_von_mises(loaded_frame, points_per_member, colormap, show_colorbar, value_limits, save_path)
    elif result_type == "deflection":
        plot_deflection(loaded_frame, scale_factor, points_per_member, colormap, show_undeformed, show_colorbar, save_path)
    else:
        raise ValueError(f"Unsupported result_type: '{result_type}'. Use 'von_mises' or 'deflection'.")


def plot_member_diagrams(
    loaded_frame: 'LoadedFrame',
    member_id: str,
    save_path: Optional[str] = None
) -> None:
    """
    Plot internal force diagrams (N, V, M, T) for a specific member.
    
    Creates a 2x2 subplot figure showing axial, shear, bending, and torsion.
    
    Args:
        loaded_frame: LoadedFrame object with analysis results
        member_id: Member identifier
        save_path: Path to save the plot
    """
    member_results = loaded_frame.get_member_results(member_id)
    member = loaded_frame.frame.get_member(member_id)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Axial force
    ax = axes[0, 0]
    axial = member_results.axial
    ax.plot([x for x, _ in axial.action], [v for _, v in axial.action], 'b-', linewidth=2)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Position along member (m)')
    ax.set_ylabel('Axial Force N (N)')
    ax.set_title('Axial Force Diagram')
    ax.grid(True, alpha=0.3)
    
    # Shear forces
    ax = axes[0, 1]
    shear_y = member_results.shear_y
    shear_z = member_results.shear_z
    ax.plot([x for x, _ in shear_y.action], [v for _, v in shear_y.action], 'r-', linewidth=2, label='Vy')
    ax.plot([x for x, _ in shear_z.action], [v for _, v in shear_z.action], 'g-', linewidth=2, label='Vz')
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Position along member (m)')
    ax.set_ylabel('Shear Force (N)')
    ax.set_title('Shear Force Diagrams')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Bending moments
    ax = axes[1, 0]
    bending_y = member_results.bending_y
    bending_z = member_results.bending_z
    ax.plot([x for x, _ in bending_y.action], [v for _, v in bending_y.action], 'r-', linewidth=2, label='My')
    ax.plot([x for x, _ in bending_z.action], [v for _, v in bending_z.action], 'g-', linewidth=2, label='Mz')
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Position along member (m)')
    ax.set_ylabel('Bending Moment (N·m)')
    ax.set_title('Bending Moment Diagrams')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Torsion
    ax = axes[1, 1]
    torsion = member_results.torsion
    ax.plot([x for x, _ in torsion.action], [v for _, v in torsion.action], 'm-', linewidth=2)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Position along member (m)')
    ax.set_ylabel('Torsion T (N·m)')
    ax.set_title('Torsion Diagram')
    ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'Internal Force Diagrams - Member {member_id}', fontsize=14, weight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
    else:
        plt.show()


def _compute_arrow_scale(loaded_frame: 'LoadedFrame') -> float:
    """Compute appropriate arrow scale based on frame size."""
    # Get bounding box of frame
    all_positions = list(loaded_frame.frame.node_positions.values())
    positions_array = np.array(all_positions)
    
    ranges = positions_array.max(axis=0) - positions_array.min(axis=0)
    max_range = np.max(ranges)
    
    # Arrow length as 10% of max dimension
    return max_range * 0.1


def _set_axes_equal(ax):
    """Set 3D plot axes to equal scale."""
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])
