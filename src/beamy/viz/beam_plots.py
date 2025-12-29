from __future__ import annotations
from typing import TYPE_CHECKING, Optional, List
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import numpy as np

if TYPE_CHECKING:
    from ..beam1d.analysis import LoadedMember
    from ..core.support import Support
    from sectiony import Section


def _normalize_svg_path(save_path: str) -> str:
    root, ext = os.path.splitext(save_path)
    if ext.lower() != ".svg":
        return f"{root}.svg"
    return save_path


def _create_arrow_cone(tip, direction, cone_length, cone_radius, num_segments=8):
    base_center = tip - direction * cone_length
    if abs(direction[0]) < 0.9: perp1 = np.array([1, 0, 0])
    else: perp1 = np.array([0, 1, 0])
    perp1 = perp1 - np.dot(perp1, direction) * direction
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(direction, perp1)
    perp2 /= np.linalg.norm(perp2)
    base_points = [base_center + cone_radius * (np.cos(2*np.pi*i/num_segments)*perp1 + np.sin(2*np.pi*i/num_segments)*perp2) for i in range(num_segments)]
    side_faces = [[tip, base_points[i], base_points[(i+1)%num_segments]] for i in range(num_segments)]
    return side_faces, [base_points]

def _plot_point_forces(ax, point_forces, beam_length):
    if not point_forces: return
    max_mag = max([np.linalg.norm(pf.force) for pf in point_forces]) or 1.0
    max_arrow_L = beam_length / 5.0
    all_s, all_b = [], []
    for pf in point_forces:
        px, py, pz = pf.point; fx, fy, fz = pf.force
        tip, f_vec = np.array([pz, px, py]), np.array([fz, fx, fy])
        mag = np.linalg.norm(f_vec)
        arrow_L = (mag / max_mag) * max_arrow_L if max_mag > 0 else max_arrow_L
        direction = f_vec / mag if mag > 0 else np.array([1,0,0])
        tail = tip - direction * arrow_L
        side, base = _create_arrow_cone(tip, direction, arrow_L*0.2, arrow_L*0.08, 12)
        all_s.extend(side); all_b.extend(base)
        c_base = tip - direction * arrow_L*0.2
        ax.plot([tail[0], c_base[0]], [tail[1], c_base[1]], [tail[2], c_base[2]], color='red', lw=2, zorder=100)
    if all_s: ax.add_collection3d(Poly3DCollection(all_s, facecolor='red', edgecolor='red'))
    if all_b: ax.add_collection3d(Poly3DCollection(all_b, facecolor='#ff4444', edgecolor='#ff4444'))

def _plot_distributed_forces(ax, dist_forces, beam_length):
    if not dist_forces: return
    max_mag = max([max(np.linalg.norm(df.start_force), np.linalg.norm(df.end_force)) for df in dist_forces]) or 1.0
    max_arrow_L = beam_length / 5.0
    all_s, all_b = [], []
    for df in dist_forces:
        L_load = np.linalg.norm(df.end_position - df.start_position)
        if L_load < 1e-10: continue
        avg_mag = (np.linalg.norm(df.start_force) + np.linalg.norm(df.end_force)) / 2.0
        avg_L = (avg_mag / max_mag) * max_arrow_L
        n_arrows = max(2, int(np.ceil((L_load / beam_length) * 10)))
        ts = np.linspace(0, 1, n_arrows)
        pos_interp = np.outer(1-ts, df.start_position) + np.outer(ts, df.end_position)
        force_interp = np.outer(1-ts, df.start_force) + np.outer(ts, df.end_force)
        tails = []
        for i in range(n_arrows):
            tip, f_vec = np.array([pos_interp[i,2], pos_interp[i,0], pos_interp[i,1]]), np.array([force_interp[i,2], force_interp[i,0], force_interp[i,1]])
            mag = np.linalg.norm(f_vec)
            if mag < 1e-10: continue
            arrow_L = (mag / max_mag) * max_arrow_L
            direction = f_vec / mag
            tail = tip - direction * arrow_L; tails.append(tail)
            side, base = _create_arrow_cone(tip, direction, avg_L*0.2, avg_L*0.08, 12)
            all_s.extend(side); all_b.extend(base)
            ax.plot([tail[0], tip[0]-direction*avg_L*0.2], [tail[1], tip[1]-direction*avg_L*0.2], [tail[2], tip[2]-direction*avg_L*0.2], color='green', lw=2, zorder=100)
        if len(tails) >= 2:
            tails = np.array(tails)
            ax.plot(tails[:,0], tails[:,1], tails[:,2], color='green', lw=2, zorder=100)
    if all_s: ax.add_collection3d(Poly3DCollection(all_s, facecolor='green', edgecolor='green'))
    if all_b: ax.add_collection3d(Poly3DCollection(all_b, facecolor='#44dd44', edgecolor='#44dd44'))

def _plot_moments(ax, moments, beam_L):
    if not moments: return
    max_mag = max([np.linalg.norm(m.moment) for m in moments]) or 1.0
    max_R = beam_L / 10.0
    all_s, all_b = [], []
    for m in moments:
        mag = np.linalg.norm(m.moment)
        if mag < 1e-10: continue
        center = np.array([0, m.x, 0])
        R = (mag / max_mag) * max_R
        m_plot = np.array([m.moment[2], m.moment[0], m.moment[1]])
        n = m_plot / np.linalg.norm(m_plot)
        u = np.array([1,0,0]) if abs(n[0]) < 0.9 else np.array([0,1,0])
        u = (u - np.dot(u, n) * n); u /= np.linalg.norm(u)
        v = np.cross(n, u)
        angles = np.linspace(0, 1.5*np.pi, 30)
        arc = np.array([center + R*(np.cos(a)*u + np.sin(a)*v) for a in angles])
        end_dir = -np.sin(1.5*np.pi)*u + np.cos(1.5*np.pi)*v
        tip = arc[-1] + end_dir * R*0.3
        side, base = _create_arrow_cone(tip, end_dir, R*0.3, R*0.12, 12)
        all_s.extend(side); all_b.extend(base)
        ax.plot(arc[:,0], arc[:,1], arc[:,2], color='blue', lw=2, zorder=100)
    if all_s: ax.add_collection3d(Poly3DCollection(all_s, facecolor='blue', edgecolor='blue'))
    if all_b: ax.add_collection3d(Poly3DCollection(all_b, facecolor='#4444dd', edgecolor='#4444dd'))

def plot_beam_diagram(loaded_member: LoadedMember, plot_stress=False, plot_section=True, plot_supports=True, save_path=None):
    # vNext: `LoadedMember` is a 2-node / 1-member frame analysis wrapper.
    # Delegate plotting to the generic frame plotter to support arbitrary member orientation.
    from .frame_plots import plot_frame, plot_von_mises

    frame = loaded_member.frame
    if plot_stress:
        plot_von_mises(frame, save_path=save_path)
        return

    plot_frame(
        frame,
        show_loads=True,
        show_moments=True,
        show_member_ids=False,
        show_node_ids=True,
        deformed=False,
        scale_factor=1.0,
        show_reactions=False,
        save_path=save_path,
    )

def plot_analysis_results(loaded_member: LoadedMember, save_path=None, show=True, points=100, units=None):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 8))
    u = units or {}
    ul, uf, um, ud = [f" ({u[k]})" if k in u else "" for k in ['length', 'force', 'moment', 'deflection']]
    sy, sz = loaded_member.shear("y", points), loaded_member.shear("z", points)
    ax1.plot(sy.action._x, sy.action._values, label="y-axis"); ax1.plot(sz.action._x, sz.action._values, label="z-axis", ls="--")
    ax1.set_title("Shear Force"); ax1.set_xlabel(f"Position{ul}"); ax1.set_ylabel(f"Force{uf}"); ax1.legend(); ax1.grid(True, ls=':', alpha=0.6); ax1.axhline(0, c='k', lw=0.5)
    by, bz = loaded_member.bending("y", points), loaded_member.bending("z", points)
    ax2.plot(by.action._x, by.action._values, label="z-axis"); ax2.plot(bz.action._x, bz.action._values, label="y-axis", ls="--")
    ax2.set_title("Bending Moment"); ax2.set_xlabel(f"Position{ul}"); ax2.set_ylabel(f"Moment{um}"); ax2.legend(); ax2.grid(True, ls=':', alpha=0.6); ax2.axhline(0, c='k', lw=0.5)
    dy, dz = loaded_member.deflection("y", points), loaded_member.deflection("z", points)
    ax3.plot(dy._x, dy._values, label="y-axis"); ax3.plot(dz._x, dz._values, label="z-axis", ls="--")
    ax3.set_title("Deflection"); ax3.set_xlabel(f"Position{ul}"); ax3.set_ylabel(f"Displacement{ud}"); ax3.legend(); ax3.grid(True, ls=':', alpha=0.6); ax3.axhline(0, c='k', lw=0.5)
    ax_res, tor = loaded_member.axial(points), loaded_member.torsion(points)
    l1 = ax4.plot(ax_res.action._x, ax_res.action._values, label="Axial Force", color="red"); ax4.set_ylabel(f"Axial Force{uf}", color="red")
    ax4_r = ax4.twinx(); l2 = ax4_r.plot(tor.action._x, tor.action._values, label="Torsion", color="purple", ls="--"); ax4_r.set_ylabel(f"Torsion Moment{um}", color="purple")
    ax4.set_title("Axial Force & Torsion"); ax4.set_xlabel(f"Position{ul}"); ax4.grid(True, ls=':', alpha=0.6); ax4.axhline(0, c='k', lw=0.5)
    lns = l1+l2; ax4.legend(lns, [l.get_label() for l in lns])
    plt.tight_layout()
    if save_path:
        save_path = _normalize_svg_path(str(save_path))
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    if show: plt.show()

class StressPlotter:
    """Plots stress distributions on the beam's cross-section."""
    def __init__(self, loaded_member: LoadedMember):
        self.loaded_member = loaded_member

    def _get_internal_forces_at(self, x_pos: float) -> dict[str, float]:
        def get_val(analysis_res, x):
            return np.interp(x, analysis_res.action._x, analysis_res.action._values)
        return {
            "N": get_val(self.loaded_member.axial(), x_pos),
            "Vy": get_val(self.loaded_member.shear("y"), x_pos),
            "Vz": get_val(self.loaded_member.shear("z"), x_pos),
            "Mx": get_val(self.loaded_member.torsion(), x_pos),
            "My": get_val(self.loaded_member.bending("z"), x_pos),
            "Mz": get_val(self.loaded_member.bending("y"), x_pos)
        }

    def plot_stress_at(self, x_pos: float, stress_type="von_mises", ax=None, show=True, cmap="viridis", title=None):
        from sectiony.stress import Stress
        f = self._get_internal_forces_at(x_pos)
        stress = Stress(
            section=self.loaded_member.section,
            N=f["N"],
            Vy=f["Vy"],
            Vz=f["Vz"],
            Mx=f["Mx"],
            My=f["My"],
            Mz=f["Mz"],
        )
        ax = stress.plot(stress_type=stress_type, ax=ax, show=False, cmap=cmap)
        if ax:
            ax.set_title(title or f"{ax.get_title()}\n@ x={x_pos:.3f}")
            if show: plt.show()
        return ax

def plot_section(section: Section, ax: Optional[plt.Axes] = None, show: bool = True) -> Optional[plt.Axes]:
    from sectiony.plotter import plot_section as sectiony_plot_section
    return sectiony_plot_section(section, ax=ax, show=show)

def plot_loads(ax, load_case, beam_length: float):
    if not load_case: return
    # Legacy helper (Beam1D-style loads) retained for backward scripts.
    if hasattr(load_case, "point_forces"):
        _plot_point_forces(ax, load_case.point_forces, beam_length)
        _plot_moments(ax, load_case.moments, beam_length)
        _plot_distributed_forces(ax, load_case.dist_forces, beam_length)

def plot_supports(supports: List[Support], beam_length: float, unit: str = "m", save_path: Optional[str] = None, show: bool = True):
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.plot([0, beam_length], [0, 0], color='black', lw=2, zorder=1)
    if supports:
        xs = [s.x for s in supports]
        ax.scatter(xs, [0]*len(xs), marker='o', color='red', s=100, zorder=2)
        for s in supports:
            ax.annotate(s.type, (s.x, 0), xytext=(0, 15), textcoords='offset points', ha='center', rotation=45)
            ax.annotate(f"{s.x} {unit}", (s.x, 0), xytext=(0, -20), textcoords='offset points', ha='center', va='top')
    ax.set_xlim(-beam_length*0.1, beam_length*1.1); ax.set_ylim(-1, 1); ax.axis('off')
    plt.tight_layout()
    if save_path:
        save_path = _normalize_svg_path(str(save_path))
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    if show: plt.show()

