from __future__ import annotations
from dataclasses import dataclass, field
from typing import Sequence, List, Tuple, Dict, Callable, Optional, Any
import numpy as np

from .beam import Beam1D
from ..core.support import Support
from ..core.results import Result, AnalysisResult
from ..core.loads import LoadCase

def _solve_fem_1d(nodes, dof_per_node, k_local_fn, f_global_fn, is_fixed_fn):
    n_nodes = len(nodes)
    ndof = dof_per_node * n_nodes
    K = np.zeros((ndof, ndof))
    x_nodes = np.array([s.x for s in nodes])
    
    for e in range(n_nodes - 1):
        L = x_nodes[e+1] - x_nodes[e]
        ke = k_local_fn(L)
        idx = []
        for i in (e, e+1):
            for k in range(dof_per_node):
                idx.append(i * dof_per_node + k)
        for a in range(len(idx)):
            for b in range(len(idx)):
                K[idx[a], idx[b]] += ke[a, b]

    f = f_global_fn(nodes, ndof)
    fixed_dofs = sorted([i*dof_per_node + k for i, s in enumerate(nodes) 
                         for k in range(dof_per_node) if is_fixed_fn(s, k)])
    free_dofs = [d for d in range(ndof) if d not in fixed_dofs]
    
    d = np.zeros(ndof)
    if free_dofs:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        f_f = f[free_dofs]
        d[free_dofs] = np.linalg.solve(K_ff, f_f)
        
    return d, K @ d - f

def solve_x_reactions(beam, loads):
    supports_sorted = sorted(beam.supports, key=lambda s: s.x)
    existing_x = {s.x for s in supports_sorted}
    
    def k_linear(L): return (1.0/L) * np.array([[1., -1.], [-1., 1.]])

    # Axial
    nodes_axial = sorted(supports_sorted + [Support(x, "000000") for x, _ in loads.Fxs if x not in existing_x], key=lambda s: s.x)
    x_nodes_axial = np.array([s.x for s in nodes_axial])
    def build_f_axial(nodes, ndof):
        f = np.zeros(ndof)
        for x_p, val in loads.Fxs:
            idx = np.searchsorted(x_nodes_axial, x_p)
            f[idx] += val
        return f
    d_x, r_x = _solve_fem_1d(nodes_axial, 1, k_linear, build_f_axial, lambda s, k: s.type[0] == "1")

    # Torsion
    t_load_x = {x for x, _ in loads.Mxs} | {x for x, _ in loads.Fys} | {x for x, _ in loads.Fzs}
    nodes_torsion = sorted(supports_sorted + [Support(x, "000000") for x in t_load_x if x not in existing_x], key=lambda s: s.x)
    x_nodes_torsion = np.array([s.x for s in nodes_torsion])
    sc_y, sc_z = getattr(beam.section, 'SCy', 0.0), getattr(beam.section, 'SCz', 0.0)
    def build_f_torsion(nodes, ndof):
        f = np.zeros(ndof)
        for x_p, val in loads.Mxs: f[np.searchsorted(x_nodes_torsion, x_p)] += val
        if sc_y or sc_z:
            for x_p, fy in loads.Fys: f[np.searchsorted(x_nodes_torsion, x_p)] += sc_z * fy
            for x_p, fz in loads.Fzs: f[np.searchsorted(x_nodes_torsion, x_p)] -= sc_y * fz
        return f
    d_rx, r_rx = _solve_fem_1d(nodes_torsion, 1, k_linear, build_f_torsion, lambda s, k: s.type[3] == "1")

    for s in supports_sorted:
        if s.x in x_nodes_axial: s.reactions["Fx"] = float(r_x[np.searchsorted(x_nodes_axial, s.x)])
        if s.x in x_nodes_torsion: s.reactions["Mx"] = float(r_rx[np.searchsorted(x_nodes_torsion, s.x)])
    return d_x, d_rx, x_nodes_axial, x_nodes_torsion

def solve_transverse_reactions(beam, loads, axis="z"):
    if axis == "z": trans_idx, rot_idx, shear_key, moment_key, I_attr, s_pairs, m_pairs = 2, 4, "Fz", "My", "Iy", loads.Fzs, loads.Mys
    else: trans_idx, rot_idx, shear_key, moment_key, I_attr, s_pairs, m_pairs = 1, 5, "Fy", "Mz", "Iz", loads.Fys, loads.Mzs
    EI = beam.material.E * getattr(beam.section, I_attr)
    supports_sorted = sorted(beam.supports, key=lambda s: s.x)
    existing_x = {s.x for s in supports_sorted}
    all_load_x = {float(x) for x, _ in s_pairs} | {float(x) for x, _ in m_pairs} | {0.0, beam.L}
    nodes = sorted(supports_sorted + [Support(x, "000000") for x in all_load_x if x not in existing_x], key=lambda s: s.x)
    x_nodes = np.array([s.x for s in nodes])
    
    def k_hermite(L):
        k = EI / L**3
        return k * np.array([[12, 6*L, -12, 6*L], [6*L, 4*L**2, -6*L, 2*L**2], [-12, -6*L, 12, -6*L], [6*L, 2*L**2, -6*L, 4*L**2]])
    def build_f(nodes, ndof):
        f = np.zeros(ndof)
        for x_p, Fw in s_pairs: f[2*np.searchsorted(x_nodes, x_p)] += Fw
        for x_p, Mth in m_pairs: f[2*np.searchsorted(x_nodes, x_p)+1] += Mth
        return f
    d, r = _solve_fem_1d(nodes, 2, k_hermite, build_f, lambda s, k: s.type[trans_idx if k==0 else rot_idx] == "1")
    for s in supports_sorted:
        idx = np.searchsorted(x_nodes, s.x)
        s.reactions[shear_key], s.reactions[moment_key] = float(r[2*idx]), float(r[2*idx+1])
    return d, x_nodes

def get_all_loads(loads, beam):
    load_map = {}
    def add(x, t, v):
        k = (float(x), t)
        if k in load_map: load_map[k] += float(v)
        else: load_map[k] = float(v)
    for x, f in loads.Fxs: add(x, "Fx", f)
    for x, f in loads.Fys: add(x, "Fy", f)
    for x, f in loads.Fzs: add(x, "Fz", f)
    for x, m in loads.Mxs: add(x, "Mx", m)
    for x, m in loads.Mys: add(x, "My", m)
    for x, m in loads.Mzs: add(x, "Mz", m)
    for s in beam.supports:
        for k, v in s.reactions.items():
            if v != 0: add(s.x, "R"+k if not k.startswith("R") else k, v)
    return sorted([(x, t, v) for (x, t), v in load_map.items()], key=lambda x: x[0])

def _accumulate_loads(points, loads, moment_loads=None):
    f_vals, m_vals = np.zeros_like(points), np.zeros_like(points)
    for i, x in enumerate(points):
        for lx, lv in loads:
            if lx <= x + 1e-9:
                f_vals[i] += lv
                if moment_loads is not None: m_vals[i] += lv * (x - lx)
        if moment_loads:
            for mx, mv in moment_loads:
                if mx <= x + 1e-9: m_vals[i] += mv
    return f_vals, m_vals

def _hermite_displacement(points, x_nodes, d_vec):
    disps = np.zeros_like(points)
    for i, x in enumerate(points):
        idx = max(0, min(len(x_nodes)-2, np.searchsorted(x_nodes, x) - 1))
        L = x_nodes[idx+1] - x_nodes[idx]
        xi = (x - x_nodes[idx]) / L
        w1, t1, w2, t2 = d_vec[2*idx:2*idx+4]
        xi2, xi3 = xi*xi, xi*xi*xi
        N1, N2, N3, N4 = 1-3*xi2+2*xi3, L*(xi-2*xi2+xi3), 3*xi2-2*xi3, L*(-xi2+xi3)
        disps[i] = N1*w1 + N2*t1 + N3*w2 + N4*t2
    return disps

@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase
    all_loads: List[Tuple[float, str, float]] = field(init=False)
    _d_x: np.ndarray = field(init=False); _d_rx: np.ndarray = field(init=False)
    _x_nodes_axial: np.ndarray = field(init=False); _x_nodes_torsion: np.ndarray = field(init=False)
    _d_y: np.ndarray = field(init=False); _x_nodes_y: np.ndarray = field(init=False)
    _d_z: np.ndarray = field(init=False); _x_nodes_z: np.ndarray = field(init=False)

    def __post_init__(self):
        self._d_x, self._d_rx, self._x_nodes_axial, self._x_nodes_torsion = solve_x_reactions(self.beam, self.loads)
        self._d_y, self._x_nodes_y = solve_transverse_reactions(self.beam, self.loads, axis="y")
        self._d_z, self._x_nodes_z = solve_transverse_reactions(self.beam, self.loads, axis="z")
        self.all_loads = get_all_loads(self.loads, self.beam)

    def shear(self, axis, points=100): return self._transverse_analysis(axis, points, "shear")
    def bending(self, axis, points=100): return self._transverse_analysis(axis, points, "bending")
    def axial(self, points=100): return self._axial_analysis(points, "axial")
    def torsion(self, points=100): return self._axial_analysis(points, "torsion")
    def deflection(self, axis, points=100): return self.bending(axis, points).displacement

    def von_mises(self, points=100):
        r_ax = self.axial(points).stress; r_by = self.bending("y", points).stress; r_bz = self.bending("z", points).stress
        sigma = np.abs(r_ax._values) + np.abs(r_by._values) + np.abs(r_bz._values)
        r_sy = self.shear("y", points).stress; r_sz = self.shear("z", points).stress; r_tor = self.torsion(points).stress
        tau = np.abs(r_sy._values) + np.abs(r_sz._values) + np.abs(r_tor._values)
        return Result(r_ax._x, np.sqrt(sigma**2 + 3 * tau**2))

    def check(self, check_module: Any, **kwargs) -> Any:
        """
        Run a check module on this loaded beam.
        The check module must have a run(loaded_beam, **kwargs) function.
        """
        if hasattr(check_module, "run"):
            return check_module.run(self, **kwargs)
        if callable(check_module):
            return check_module(self, **kwargs)
        raise ValueError("check_module must be a module with a run() function or a callable.")

    def check_aisc_chapter_f(self, length_unit, force_unit) -> Any:
        from ..checks import aisc_9
        return self.check(aisc_9, length_unit=length_unit, force_unit=force_unit)

    def plot(self, **kwargs):
        from ..viz.beam_plots import plot_beam_diagram
        plot_beam_diagram(self, **kwargs)

    def plot_results(self, **kwargs):
        from ..viz.beam_plots import plot_analysis_results
        plot_analysis_results(self, **kwargs)

    def _transverse_analysis(self, axis, points, mode):
        xs = np.linspace(0, self.beam.L, points)
        if axis == "y": s_k, m_k, d_v, x_n, I, c = "Fy", "Mz", self._d_y, self._x_nodes_y, self.beam.section.Iz, self.beam.section.y_max
        else: s_k, m_k, d_v, x_n, I, c = "Fz", "My", self._d_z, self._x_nodes_z, self.beam.section.Iy, self.beam.section.z_max
        A = self.beam.section.A
        r_f = [(x, v) for x, t, v in self.all_loads if t in (s_k, "R"+s_k[-1])]
        r_m = [(x, v) for x, t, v in self.all_loads if t in (m_k, "RM"+m_k[-1])]
        V, M = _accumulate_loads(xs, r_f, r_m)
        disps = _hermite_displacement(xs, x_n, d_v) if d_v is not None else np.zeros(points)
        if mode == "shear": return AnalysisResult(Result(xs, V), Result(xs, V/A if A>0 else V*0), Result(xs, disps))
        return AnalysisResult(Result(xs, M), Result(xs, M*c/I if I>0 else M*0), Result(xs, disps))

    def _axial_analysis(self, points, mode):
        xs = np.linspace(0, self.beam.L, points)
        if mode == "axial": k, r_k, prop, x_n, d_v, f = "Fx", "Rx", self.beam.section.A, self._x_nodes_axial, self._d_x, lambda p: 1.0/p if p>0 else 0.0
        else: k, r_k, prop, x_n, d_v, f = "Mx", "RMx", self.beam.section.J, self._x_nodes_torsion, self._d_rx, lambda p: max(abs(self.beam.section.y_max), abs(self.beam.section.z_max))/p if p>0 else 0.0
        rel = [(x, v) for x, t, v in self.all_loads if t in (k, r_k)]
        actions, _ = _accumulate_loads(xs, rel)
        disps = np.interp(xs, x_n, d_v) if d_v is not None else np.zeros_like(xs)
        return AnalysisResult(Result(xs, actions), Result(xs, actions * f(prop)), Result(xs, disps))
