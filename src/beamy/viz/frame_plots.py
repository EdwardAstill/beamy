from __future__ import annotations
from typing import Optional, Tuple, TYPE_CHECKING
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

if TYPE_CHECKING:
    from ..frame.analysis import LoadedFrame
    from ..frame.member import Member


def _beam_shape_functions(xi: float, L: float) -> tuple[float, float, float, float]:
    """
    Cubic Euler-Bernoulli beam shape functions for transverse deflection.

    Returns (N1, N2, N3, N4) such that:
      w(x) = N1*w1 + N2*theta1 + N3*w2 + N4*theta2

    where w is transverse displacement and theta is the corresponding rotation DOF.
    """
    n1 = 1.0 - 3.0 * xi**2 + 2.0 * xi**3
    n2 = L * (xi - 2.0 * xi**2 + xi**3)
    n3 = 3.0 * xi**2 - 2.0 * xi**3
    n4 = L * (-xi**2 + xi**3)
    return n1, n2, n3, n4


def _member_deformed_points(
    loaded_frame: LoadedFrame,
    member: "Member",
    points_per_member: int,
    scale_factor: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute smooth (curved) deformed coordinates along a member using
    Euler-Bernoulli shape functions (uses translations + rotations).

    Returns (pts_global, disp_mags) where:
      - pts_global: (n,3) array of deformed global points
      - disp_mags:  (n,) array of displacement magnitudes at each point (global)
    """
    node_pos = loaded_frame.frame.node_positions
    s = node_pos[member.start_node_id]
    e = node_pos[member.end_node_id]
    L = float(np.linalg.norm(e - s))
    if L < 1e-12:
        pts = np.repeat(s.reshape(1, 3), repeats=points_per_member, axis=0)
        mags = np.zeros(points_per_member)
        return pts, mags

    # Global nodal DOFs: [UX, UY, UZ, RX, RY, RZ]
    d_s_g = loaded_frame.nodal_displacements[member.start_node_id]
    d_e_g = loaded_frame.nodal_displacements[member.end_node_id]

    # Truss/cable plotting:
    # These elements should NOT inherit beam end-rotations in the visualization.
    # Draw them as straight lines between displaced nodes (translations only).
    if getattr(member, "element_type", "beam") != "beam":
        xis = np.linspace(0.0, 1.0, points_per_member)
        pts_global = np.zeros((points_per_member, 3))
        disp_mags = np.zeros(points_per_member)

        u_s = d_s_g[:3]
        u_e = d_e_g[:3]
        for i, xi in enumerate(xis):
            p0 = (1.0 - xi) * s + xi * e
            u = (1.0 - xi) * u_s + xi * u_e
            pts_global[i, :] = p0 + scale_factor * u
            disp_mags[i] = float(np.linalg.norm(u))

        return pts_global, disp_mags

    # Transform to member local coordinates
    # member.transformation_matrix has rows = local axes in global components
    r = member.transformation_matrix
    u_s_l = r @ d_s_g[:3]
    u_e_l = r @ d_e_g[:3]
    rot_s_l = r @ d_s_g[3:]
    rot_e_l = r @ d_e_g[3:]

    # Local DOFs
    u1 = float(u_s_l[0])
    v1 = float(u_s_l[1])
    w1 = float(u_s_l[2])
    rx1 = float(rot_s_l[0])
    ry1 = float(rot_s_l[1])
    rz1 = float(rot_s_l[2])

    u2 = float(u_e_l[0])
    v2 = float(u_e_l[1])
    w2 = float(u_e_l[2])
    rx2 = float(rot_e_l[0])
    ry2 = float(rot_e_l[1])
    rz2 = float(rot_e_l[2])

    xis = np.linspace(0.0, 1.0, points_per_member)
    pts_global = np.zeros((points_per_member, 3))
    disp_mags = np.zeros(points_per_member)

    for i, xi in enumerate(xis):
        # Undeformed global point
        p0 = (1.0 - xi) * s + xi * e

        # Axial + torsion (linear interpolation)
        u = (1.0 - xi) * u1 + xi * u2
        # torsion phi not used for plotting centerline, but keep for completeness
        _phi = (1.0 - xi) * rx1 + xi * rx2

        # Transverse bending (cubic)
        n1, n2, n3, n4 = _beam_shape_functions(xi, L)
        # v (local y displacement) uses rotation about local z (rz).
        # A positive rz (right-hand rule about z) produces a positive dv/dx.
        v = n1 * v1 + n2 * rz1 + n3 * v2 + n4 * rz2
        
        # w (local z displacement) uses rotation about local y (ry).
        # In the local x-z plane, a positive ry (right-hand rule about y) 
        # actually produces a NEGATIVE dw/dx slope (moving from Z towards X).
        w = n1 * w1 + n2 * (-ry1) + n3 * w2 + n4 * (-ry2)

        d_l = np.array([u, v, w], dtype=float)
        d_g = r.T @ d_l

        pts_global[i, :] = p0 + scale_factor * d_g
        disp_mags[i] = float(np.linalg.norm(d_g))

    return pts_global, disp_mags

def _create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=8):
    """Create 3D arrow cone for force visualization."""
    base_center = tip - direction * cone_length
    perp1 = np.array([1, 0, 0]) if abs(direction[0]) < 0.9 else np.array([0, 1, 0])
    perp1 = perp1 - np.dot(perp1, direction) * direction
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(direction, perp1)
    perp2 /= np.linalg.norm(perp2)
    base_pts = [base_center + cone_radius * (np.cos(2*np.pi*i/num_segments)*perp1 + np.sin(2*np.pi*i/num_segments)*perp2) for i in range(num_segments)]
    side_faces = [[tip, base_pts[i], base_pts[(i+1)%num_segments]] for i in range(num_segments)]
    return side_faces, [base_pts]

def plot_frame(
    loaded_frame: LoadedFrame,
    show_loads: bool = True,
    show_moments: bool = True,
    show_member_ids: bool = False,
    show_node_ids: bool = True,
    deformed: bool = False,
    scale_factor: float = 1.0,
    show_reactions: bool = False,
    save_path: Optional[str] = None
) -> None:
    """Simplified frame geometry plot with constraint labels, force arrows, and moment vectors."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    node_pos = loaded_frame.frame.node_positions
    
    # Plot members (red if member constraints are present, black otherwise)
    for member in loaded_frame.frame.members:
        s, e = node_pos[member.start_node_id], node_pos[member.end_node_id]
        is_constrained = bool(member.constraints) and ("1" in member.constraints)
        color = 'red' if is_constrained else 'black'
        ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], color=color, lw=2)
    
    # Plot deformed shape if requested (use smooth shape functions)
    if deformed:
        for member in loaded_frame.frame.members:
            pts, _mags = _member_deformed_points(
                loaded_frame=loaded_frame,
                member=member,
                points_per_member=25,
                scale_factor=scale_factor,
            )
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], "b--", lw=1.5, alpha=0.6, zorder=5)
    
    # Plot nodes and constraints
    for nid, pos in node_pos.items():
        node = loaded_frame.frame.nodes[nid]
        if node.support:
            # Constrained node: red dot with constraint label (inline)
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='red', marker='o', s=80, zorder=10)
            # Format constraint as 6-digit string (translations + rotations)
            constraint_label = str(node.support)
            ax.text(pos[0], pos[1], pos[2], f'  {constraint_label}', fontsize=9, color='red', weight='bold')
        else:
            # Free node: small black dot
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='black', marker='o', s=40, zorder=10)
            
            # Label free nodes only (constrained nodes already have their support code displayed)
            if show_node_ids:
                offset = pos[2] * 0.08
                ax.text(pos[0], pos[1], pos[2] + offset, nid, fontsize=8, ha='center', color='gray')
    
    # Plot member IDs at midpoints
    if show_member_ids:
        for m in loaded_frame.frame.members:
            mid = (node_pos[m.start_node_id] + node_pos[m.end_node_id]) / 2
            ax.text(mid[0], mid[1], mid[2], m.id, fontsize=8, color='blue', style='italic', alpha=0.7)
    
    # Plot applied forces with arrow cones (like 1D beam analysis)
    if show_loads:
        all_forces = [(node_pos[nf.node_id], nf.force) for nf in loaded_frame.loads.nodal_forces]
        if all_forces:
            max_force = max([np.linalg.norm(f) for _, f in all_forces])
            if max_force > 1e-10:
                # Compute arrow scale based on frame size
                all_pos = np.array(list(node_pos.values()))
                frame_size = np.max(all_pos.max(axis=0) - all_pos.min(axis=0))
                arrow_scale = frame_size * 0.2
                
                all_sides, all_bases = [], []
                for pos, force in all_forces:
                    mag = np.linalg.norm(force)
                    if mag < 1e-10: continue
                    direction = force / mag
                    arrow_len = (mag / max_force) * arrow_scale
                    tip = pos
                    tail = tip - direction * arrow_len
                    
                    # Draw arrow shaft
                    cone_base = tip - direction * arrow_len * 0.25
                    ax.plot([tail[0], cone_base[0]], [tail[1], cone_base[1]], [tail[2], cone_base[2]], 
                           color='red', lw=2.5, zorder=100)
                    
                    # Draw arrow cone
                    sides, bases = _create_arrow_cone(tip, direction, arrow_len*0.25, arrow_len*0.1, 10)
                    all_sides.extend(sides)
                    all_bases.extend(bases)
                
                if all_sides:
                    ax.add_collection3d(Poly3DCollection(all_sides, facecolor='red', edgecolor='red', alpha=0.9, zorder=100))
                if all_bases:
                    ax.add_collection3d(Poly3DCollection(all_bases, facecolor='darkred', edgecolor='darkred', alpha=0.9, zorder=100))

        # Plot applied nodal moments as curved arc arrows (like 1D beam plots)
        if show_moments:
            all_moments = [(node_pos[nm.node_id], nm.moment) for nm in loaded_frame.loads.nodal_moments]
            if all_moments:
                max_moment = max([np.linalg.norm(m) for _, m in all_moments])
                if max_moment > 1e-12:
                    all_pos = np.array(list(node_pos.values()))
                    frame_size = np.max(all_pos.max(axis=0) - all_pos.min(axis=0))
                    max_R = frame_size * 0.10  # radius of arc

                    all_sides_m, all_bases_m = [], []
                    for center, moment in all_moments:
                        mag = np.linalg.norm(moment)
                        if mag < 1e-12:
                            continue
                        
                        # Arc radius scales with moment magnitude
                        R = (mag / max_moment) * max_R
                        
                        # Moment axis (right-hand rule direction)
                        n = moment / mag
                        
                        # Create two perpendicular vectors in the plane perpendicular to n
                        if abs(n[0]) < 0.9:
                            u = np.array([1.0, 0.0, 0.0])
                        else:
                            u = np.array([0.0, 1.0, 0.0])
                        u = u - np.dot(u, n) * n
                        u /= np.linalg.norm(u)
                        v = np.cross(n, u)
                        
                        # Create arc (270 degrees, 1.5*pi radians)
                        angles = np.linspace(0, 1.5 * np.pi, 30)
                        arc = np.array([center + R * (np.cos(a) * u + np.sin(a) * v) for a in angles])
                        
                        # Tangent direction at end of arc
                        end_tangent = -np.sin(1.5 * np.pi) * u + np.cos(1.5 * np.pi) * v
                        
                        # Arrow tip extends beyond arc end
                        tip = arc[-1] + end_tangent * R * 0.3
                        
                        # Draw arc
                        ax.plot(arc[:, 0], arc[:, 1], arc[:, 2], color='purple', lw=2.5, zorder=101)
                        
                        # Arrow cone at tip
                        sides, bases = _create_arrow_cone(tip, end_tangent, R * 0.3, R * 0.12, 10)
                        all_sides_m.extend(sides)
                        all_bases_m.extend(bases)

                    if all_sides_m:
                        ax.add_collection3d(
                            Poly3DCollection(all_sides_m, facecolor='purple', edgecolor='purple', alpha=0.85, zorder=101)
                        )
                    if all_bases_m:
                        ax.add_collection3d(
                            Poly3DCollection(all_bases_m, facecolor='indigo', edgecolor='indigo', alpha=0.85, zorder=101)
                        )
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Frame Geometry')
    _set_axes_equal(ax)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
    else:
        plt.show()

def plot_deflection(
    loaded_frame: LoadedFrame, scale_factor: float = 1.0, points_per_member: int = 20,
    colormap: str = "viridis", show_undeformed: bool = True, show_colorbar: bool = True, save_path: Optional[str] = None
) -> None:
    """Plot deformed frame colored by displacement magnitude."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    node_pos = loaded_frame.frame.node_positions
    
    # Collect displacement magnitudes for color mapping (use smooth interpolation)
    all_disps: list[float] = []
    for m in loaded_frame.frame.members:
        _pts, mags = _member_deformed_points(
            loaded_frame=loaded_frame,
            member=m,
            points_per_member=points_per_member,
            scale_factor=1.0,  # magnitudes should be unscaled physical displacements
        )
        all_disps.extend([float(x) for x in mags])
    
    norm = Normalize(vmin=min(all_disps), vmax=max(all_disps))
    cmap = plt.get_cmap(colormap)
    
    # Plot undeformed frame
    if show_undeformed:
        for m in loaded_frame.frame.members:
            s, e = node_pos[m.start_node_id], node_pos[m.end_node_id]
            ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], 'k--', lw=1, alpha=0.3)
    
    # Plot deformed frame with color mapping (smooth member curves)
    for m in loaded_frame.frame.members:
        pts, mags = _member_deformed_points(
            loaded_frame=loaded_frame,
            member=m,
            points_per_member=points_per_member,
            scale_factor=scale_factor,
        )
        colors = [cmap(norm(float(mag))) for mag in mags]
        for i in range(len(pts) - 1):
            ax.plot(
                [pts[i, 0], pts[i + 1, 0]],
                [pts[i, 1], pts[i + 1, 1]],
                [pts[i, 2], pts[i + 1, 2]],
                color=colors[i],
                lw=2.5,
            )
    
    if show_colorbar:
        cbar = plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), ax=ax, shrink=0.5, pad=0.1)
        cbar.set_label('Displacement Magnitude', rotation=270, labelpad=15)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Deflection (scale={scale_factor}x)')
    _set_axes_equal(ax)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
    else:
        plt.show()

def plot_von_mises(
    loaded_frame: LoadedFrame, points_per_member: int = 50, colormap: str = "turbo",
    show_colorbar: bool = True, stress_limits: Optional[Tuple[float, float]] = None, save_path: Optional[str] = None
) -> None:
    """Plot frame colored by Von Mises stress.
    
    Uses original (unsplit) members to ensure continuous stress distribution
    even when members are split at intermediate nodes.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    node_pos = loaded_frame.frame.node_positions

    # Build a single 3D line collection for performance and to avoid slow/fragile autoscaling.
    segments_list = []
    seg_stresses_list = []

    # Iterate over ORIGINAL members (not segments) for continuous stress plots
    for m in loaded_frame.original_frame.members:
        # Get member start/end positions from original frame
        s = np.asarray(loaded_frame.original_frame.get_node(m.start_node_id).position, dtype=float)
        e = np.asarray(loaded_frame.original_frame.get_node(m.end_node_id).position, dtype=float)

        # Get continuous action profile using demand_provider which stitches segments
        try:
            profile = loaded_frame.demand_provider.actions(m.id, points=points_per_member)
            xs = profile.axial._x
            
            # Compute von Mises stress from the profile actions
            # Section properties
            A = m.section.A
            Iy = m.section.Iy
            Iz = m.section.Iz
            J = m.section.J
            y_max = m.section.y_max
            z_max = m.section.z_max
            r_max = max(abs(y_max), abs(z_max))
            
            N = profile.axial._values
            Vy = profile.shear_y._values
            Vz = profile.shear_z._values
            T = profile.torsion._values
            My = profile.bending_y._values
            Mz = profile.bending_z._values
            
            # Axial stress
            sigma_axial = N / A if A > 0 else np.zeros_like(N)
            # Bending stress (max at extreme fiber)
            sigma_bending_y = np.abs(My) * z_max / Iy if Iy > 0 else np.zeros_like(My)
            sigma_bending_z = np.abs(Mz) * y_max / Iz if Iz > 0 else np.zeros_like(Mz)
            sigma = np.abs(sigma_axial) + sigma_bending_y + sigma_bending_z
            # Shear stress
            tau_shear = (np.abs(Vy) + np.abs(Vz)) / A if A > 0 else np.zeros_like(Vy)
            tau_torsion = np.abs(T) * r_max / J if J > 0 else np.zeros_like(T)
            tau = tau_shear + tau_torsion
            # Von Mises
            stresses = np.sqrt(sigma**2 + 3 * tau**2)
        except Exception as e:
            # Skip members that fail (e.g., cables, trusses)
            continue

        t_vals = xs / float(m.length) if float(m.length) > 0 else np.linspace(0.0, 1.0, len(xs))
        pts = (1.0 - t_vals)[:, None] * s + t_vals[:, None] * e

        # (n-1, 2, 3) segments
        segments_list.append(np.stack([pts[:-1], pts[1:]], axis=1))
        seg_stresses_list.append(0.5 * (stresses[:-1] + stresses[1:]))

    if not segments_list:
        return

    segments = np.concatenate(segments_list, axis=0)
    seg_stresses = np.concatenate(seg_stresses_list, axis=0)
    finite_stress_mask = np.isfinite(seg_stresses)

    if stress_limits is not None:
        vmin, vmax = float(stress_limits[0]), float(stress_limits[1])
    else:
        finite_vals = seg_stresses[finite_stress_mask]
        if finite_vals.size:
            vmin, vmax = float(np.min(finite_vals)), float(np.max(finite_vals))
        else:
            vmin, vmax = 0.0, 1.0

    if not np.isfinite(vmin):
        vmin = 0.0
    if (not np.isfinite(vmax)) or vmax <= vmin:
        vmax = vmin + 1.0

    norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap = plt.get_cmap(colormap)

    safe_stress = np.where(finite_stress_mask, seg_stresses, vmin)
    colors = cmap(norm(safe_stress))

    lc = Line3DCollection(segments, colors=colors, linewidths=3.5)
    ax.add_collection3d(lc)

    # Set axis limits from node coordinates (prevents inf limits from bad autoscale).
    pts_all = np.array(list(node_pos.values()), dtype=float)
    finite_pts = pts_all[np.all(np.isfinite(pts_all), axis=1)]
    if finite_pts.size:
        mins = finite_pts.min(axis=0)
        maxs = finite_pts.max(axis=0)
        span = maxs - mins
        pad = float(0.05 * np.max(span)) if float(np.max(span)) > 0 else 1.0
        ax.set_xlim(float(mins[0] - pad), float(maxs[0] + pad))
        ax.set_ylim(float(mins[1] - pad), float(maxs[1] + pad))
        ax.set_zlim(float(mins[2] - pad), float(maxs[2] + pad))
    
    if show_colorbar:
        cbar = plt.colorbar(ScalarMappable(cmap=cmap, norm=norm), ax=ax, shrink=0.5, pad=0.1)
        cbar.set_label('Von Mises Stress (Pa)', rotation=270, labelpad=15)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Von Mises Stress Distribution')
    _set_axes_equal(ax)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
    else:
        plt.show()

def plot_results(loaded_frame: LoadedFrame, result_type: str = "von_mises", **kwargs) -> None:
    """Unified results plot."""
    if result_type == "von_mises":
        plot_von_mises(loaded_frame, **kwargs)
    elif result_type == "deflection":
        plot_deflection(loaded_frame, **kwargs)
    else:
        raise ValueError(f"Invalid result_type: {result_type}")

def plot_member_diagrams(loaded_frame: LoadedFrame, member_id: str, save_path: Optional[str] = None) -> None:
    """Plot internal force diagrams for a specific member using direct frame analysis."""
    # Use demand_provider for continuous results across split members
    profile = loaded_frame.demand_provider.actions(member_id, points=201)
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    
    def _plot(ax, xs, vals, title, ylabel, label=None, c='b'):
        ax.plot(xs, vals, c, lw=2, label=label)
        ax.axhline(0, color='k', ls='--', lw=0.5, alpha=0.5)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Position (m)')
        ax.grid(True, alpha=0.3)
    
    _plot(axes[0,0], profile.axial._x, profile.axial._values, 'Axial Force', 'N (N)')
    _plot(axes[0,1], profile.shear_y._x, profile.shear_y._values, 'Shear Force', 'V (N)', 'Vy', 'r')
    _plot(axes[0,1], profile.shear_z._x, profile.shear_z._values, '', '', 'Vz', 'g')
    axes[0,1].legend()
    _plot(axes[1,0], profile.bending_y._x, profile.bending_y._values, 'Bending Moment', 'M (N·m)', 'My', 'r')
    _plot(axes[1,0], profile.bending_z._x, profile.bending_z._values, '', '', 'Mz', 'g')
    axes[1,0].legend()
    _plot(axes[1,1], profile.torsion._x, profile.torsion._values, 'Torsion', 'T (N·m)', c='m')
    
    fig.suptitle(f'Member {member_id} - Internal Forces', fontsize=13, weight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
    else:
        plt.show()


def plot_aisc_utilization(
    loaded_frame: LoadedFrame,
    show_node_ids: bool = True,
    save_path: Optional[str] = None,
    colormap: str = "RdYlGn_r",
    **kwargs,  # Accept but ignore legacy options like utilization_source
) -> None:
    """
    Plot frame with members colored by AISC utilization ratio.
    
    Uses direct frame analysis (MemberActionProfile) for consistent and correct results.
    Green = low utilization, Yellow = moderate, Red = high (approaching 1.0).
    For split members, uses the utilization of the original member (combined segments).
    Only one label is shown per original member (at the midpoint of the full member).
    """
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # Always use direct frame analysis method (get_aisc_utilizations)
    # This uses MemberActionProfile from demand_provider, NOT to_loaded_beam() re-solve
    utilizations = loaded_frame.get_aisc_utilizations()
    
    if not utilizations:
        print("No utilizations computed")
        return
    
    # Normalize utilizations for colormap
    util_values = [v for v in utilizations.values() if v > 0]
    max_util = max(util_values) if util_values else 1.0
    
    # Use max(max_util, 1.0) to ensure the scale goes at least to 1.0
    norm_max = max(max_util, 1.0)
    norm = Normalize(vmin=0.0, vmax=norm_max)
    cmap = plt.get_cmap(colormap)
    
    node_pos = loaded_frame.frame.node_positions
    
    # Track which original members have been labeled
    labeled_parents = set()
    
    # Plot members with color based on utilization
    for member in loaded_frame.frame.members:
        seg_id = member.id
        s, e = node_pos[member.start_node_id], node_pos[member.end_node_id]
        
        # Determine the parent member ID (for bundled/split members)
        parent_id = loaded_frame._member_segment_parent.get(seg_id, seg_id)
        
        # Look up utilization from the utilizations dict
        util = utilizations.get(parent_id, 0.0)
        
        color = cmap(norm(util))
        
        ax.plot([s[0], e[0]], [s[1], e[1]], [s[2], e[2]], color=color, lw=3, alpha=0.8)
        
        # Add text label only once per original member (at full member midpoint)
        if parent_id not in labeled_parents:
            labeled_parents.add(parent_id)
            # Get full member midpoint using original frame if available
            try:
                parent = loaded_frame.original_frame.get_member(parent_id)
                parent_start = loaded_frame.original_frame.get_node(parent.start_node_id).position
                parent_end = loaded_frame.original_frame.get_node(parent.end_node_id).position
                mid = (parent_start + parent_end) / 2
            except:
                mid = (s + e) / 2
            ax.text(mid[0], mid[1], mid[2], f'{util:.3f}', fontsize=8, ha='center',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
    # Plot nodes and constraints
    for nid, pos in node_pos.items():
        node = loaded_frame.frame.nodes[nid]
        if node.support:
            # Constrained node: red dot
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='red', marker='s', s=100, zorder=10)
        else:
            # Free node: small black dot
            ax.scatter([pos[0]], [pos[1]], [pos[2]], c='black', marker='o', s=50, zorder=10)
        
        if show_node_ids:
            offset = pos[2] * 0.05 if node.support else pos[2] * 0.08
            label = str(node.support) if node.support else nid
            ax.text(pos[0], pos[1], pos[2] + offset, label, fontsize=8, ha='center', color='gray')
    
    # Add colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('AISC Utilization Ratio', fontsize=11, weight='bold')
    
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title('AISC Chapter F Utilization Plot', fontsize=14, weight='bold')
    
    _set_axes_equal(ax)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved AISC utilization plot to {save_path}")
    else:
        plt.show()

def _set_axes_equal(ax):
    """Set 3D axes to equal scale."""
    lims = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    origin = np.mean(lims, axis=1)
    radius = 0.5 * np.max(np.abs(lims[:, 1] - lims[:, 0]))
    ax.set_xlim3d([origin[0]-radius, origin[0]+radius])
    ax.set_ylim3d([origin[1]-radius, origin[1]+radius])
    ax.set_zlim3d([origin[2]-radius, origin[2]+radius])
