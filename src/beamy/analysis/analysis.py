# analysis.py

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Sequence, List, Tuple, Dict, Callable, Optional
import numpy as np

from ..setup import Beam1D, Support
from ..setup import LoadCase


# -------------------------------------------------
# Generic 1D FEM Solver
# -------------------------------------------------
def _solve_fem_1d(
    nodes: List[Support],
    dof_per_node: int,
    k_local_fn: Callable[[float], np.ndarray],
    f_global_fn: Callable[[List[Support], int], np.ndarray],
    is_fixed_fn: Callable[[Support, int], bool],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generic solver for 1D beam problems (Axial, Torsion, Bending).

    Args:
        nodes: Sorted list of Support objects.
        dof_per_node: Number of DOFs per node (1 or 2).
        k_local_fn: Function(L) -> returns local stiffness matrix (dof_per_node*2)^2.
        f_global_fn: Function(nodes, ndof) -> returns global load vector f.
        is_fixed_fn: Function(support, local_dof_index) -> returns True if fixed.

    Returns:
        d: Global displacement vector (ndof).
        r: Global reaction vector (ndof) [r = K*d - f].
    """
    n_nodes = len(nodes)
    ndof = dof_per_node * n_nodes
    
    # 1. Build Global Stiffness Matrix K
    K = np.zeros((ndof, ndof), dtype=float)
    x_nodes = np.array([s.x for s in nodes], dtype=float)
    
    for e in range(n_nodes - 1):
        x1 = x_nodes[e]
        x2 = x_nodes[e + 1]
        L = x2 - x1
        if L <= 0.0:
            raise ValueError("Support positions must be strictly increasing in x.")
            
        ke = k_local_fn(L)
        
        # Assembly indices
        # Element connects node e and e+1
        # Global DOFs: [node_e_dofs, node_e+1_dofs]
        idx = []
        for node_idx in (e, e + 1):
            start = node_idx * dof_per_node
            for k in range(dof_per_node):
                idx.append(start + k)
        
        # Add to global K
        dim = len(idx)
        for a in range(dim):
            for b in range(dim):
                K[idx[a], idx[b]] += ke[a, b]

    # 2. Build Global Load Vector f
    f = f_global_fn(nodes, ndof)
    
    # 3. Apply Boundary Conditions
    fixed_dofs = []
    for i, s in enumerate(nodes):
        for k in range(dof_per_node):
            if is_fixed_fn(s, k):
                fixed_dofs.append(i * dof_per_node + k)
                
    fixed_dofs = np.array(sorted(list(set(fixed_dofs))), dtype=int)
    all_dofs = np.arange(ndof, dtype=int)
    free_dofs = np.array([d for d in all_dofs if d not in fixed_dofs], dtype=int)
    
    # 4. Solve for Displacements d
    d = np.zeros(ndof, dtype=float)
    
    if free_dofs.size > 0:
        K_ff = K[np.ix_(free_dofs, free_dofs)]
        f_f = f[free_dofs]
        d_f = np.linalg.solve(K_ff, f_f)
        d[free_dofs] = d_f
        
    # 5. Compute Reactions r = K*d - f
    q = K @ d
    r = q - f
    
    return d, r


# -------------------------------------------------
# Axial + torsional reactions (Fx + Mx)
# -------------------------------------------------
def solve_x_reactions(
    supports: Sequence[Support],
    loads: LoadCase,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute reactions in the x-direction DOFs (Axial Fx and Torsion Mx).
    Updates supports in place and returns global displacement vectors.
    
    Returns:
        d_x: Global axial displacements (size n_supports)
        d_rx: Global torsional displacements (size n_supports)
    """
    if not supports:
        raise ValueError("Beam has no axial/torsional supports. Unstable.")

    # Reset existing reactions
    for s in supports:
        s.reactions["Fx"] = 0.0
        s.reactions["Mx"] = 0.0

    supports_sorted = sorted(supports, key=lambda s: s.x)
    n_nodes = len(supports_sorted)
    x_nodes = np.array([s.x for s in supports_sorted], dtype=float)

    # --- Helpers for _solve_fem_1d ---
    
    # 1. Stiffness (EA or GJ factored out, so just 1/L)
    def k_linear(L: float) -> np.ndarray:
        k = 1.0 / L
        return k * np.array([[1.0, -1.0], [-1.0, 1.0]])
    
    # 2. Load Builder (Linear interpolation to nodes)
    def build_load_fn(pairs: List[Tuple[float, float]]):
        def f_fn(nodes, ndof):
            f = np.zeros(ndof, dtype=float)
            for x_p, val in pairs:
                idx = np.searchsorted(x_nodes, x_p) - 1
                if idx < 0:
                    f[0] += val
                elif idx >= n_nodes - 1:
                    f[-1] += val
                else:
                    xL, xR = x_nodes[idx], x_nodes[idx+1]
                    t = (x_p - xL) / (xR - xL)
                    f[idx] += val * (1 - t)
                    f[idx+1] += val * t
            return f
        return f_fn

    # 3. BC Checkers
    # Axial (Ux is index 0)
    def is_fixed_axial(s: Support, k: int) -> bool:
        return s.type[0] == "1"
    
    # Torsion (Rx is index 3)
    def is_fixed_torsion(s: Support, k: int) -> bool:
        return s.type[3] == "1"

    # --- Solve Axial (Fx) ---
    d_x, r_x = _solve_fem_1d(
        supports_sorted, 
        dof_per_node=1, 
        k_local_fn=k_linear,
        f_global_fn=build_load_fn(loads.Fxs),
        is_fixed_fn=is_fixed_axial
    )
    
    # --- Solve Torsion (Mx) ---
    d_rx, r_rx = _solve_fem_1d(
        supports_sorted, 
        dof_per_node=1, 
        k_local_fn=k_linear,
        f_global_fn=build_load_fn(loads.Mxs),
        is_fixed_fn=is_fixed_torsion
    )

    # Write reactions back
    for i, s in enumerate(supports_sorted):
        s.reactions["Fx"] = float(r_x[i])
        s.reactions["Mx"] = float(r_rx[i])
        
    return d_x, d_rx


# -------------------------------------------------
# Transverse reactions (Fy/Fz + My/Mz)
# -------------------------------------------------
def solve_transverse_reactions(beam: Beam1D, loads: LoadCase, axis: str = "z") -> np.ndarray:
    """
    Solve for transverse reactions using Euler–Bernoulli beam theory.
    """
    if axis not in ("y", "z"):
        raise ValueError(f"axis must be 'y' or 'z', got {axis!r}")

    # Setup parameters
    if axis == "z":
        # Bending in x-z plane (local y) -> Fz, My
        # DOFs: Uz (idx 2), Ry (idx 4)
        trans_idx, rot_idx = 2, 4
        shear_key, moment_key = "Fz", "My"
        I_attr = "Iy"
        shear_pairs, moment_pairs = loads.Fzs, loads.Mys
    else:
        # Bending in x-y plane (local z) -> Fy, Mz
        # DOFs: Uy (idx 1), Rz (idx 5)
        trans_idx, rot_idx = 1, 5
        shear_key, moment_key = "Fy", "Mz"
        I_attr = "Iz"
        shear_pairs, moment_pairs = loads.Fys, loads.Mzs

    E = beam.material.E
    I = getattr(beam.section, I_attr, None)
    if I is None:
        raise AttributeError(f"Section missing {I_attr}")
    EI = E * I

    if not beam.supports:
        raise ValueError("Beam has no supports.")

    supports_sorted = sorted(beam.supports, key=lambda s: s.x)
    x_nodes = np.array([s.x for s in supports_sorted], dtype=float)
    n_nodes = len(supports_sorted)
    dof_per_node = 2

    # Reset reactions
    for s in supports_sorted:
        s.reactions[shear_key] = 0.0
        s.reactions[moment_key] = 0.0

    # --- Helpers for _solve_fem_1d ---

    # 1. Stiffness (Hermite Cubic)
    def k_hermite(L: float) -> np.ndarray:
        L2 = L * L
        L3 = L2 * L
        k = EI / L3
        return k * np.array([
            [ 12.0,     6.0*L,   -12.0,     6.0*L],
            [  6.0*L,   4.0*L2,   -6.0*L,   2.0*L2],
            [-12.0,    -6.0*L,    12.0,    -6.0*L],
            [  6.0*L,   2.0*L2,   -6.0*L,   4.0*L2],
        ])

    # 2. Load Builder (Hermite interpolation)
    def build_load_vector(nodes, ndof) -> np.ndarray:
        f = np.zeros(ndof, dtype=float)
        
        def add_load(x_p, Fw, Mtheta):
            idx = np.searchsorted(x_nodes, x_p) - 1
            if idx < 0:
                node = 0
                f[2*node] += Fw
                f[2*node+1] += Mtheta
            elif idx >= n_nodes - 1:
                node = n_nodes - 1
                f[2*node] += Fw
                f[2*node+1] += Mtheta
            else:
                xL, xR = x_nodes[idx], x_nodes[idx+1]
                t = (x_p - xL) / (xR - xL)
                wL, wR = 1.0 - t, t
                # Force -> shear
                f[2*idx] += wL * Fw
                f[2*(idx+1)] += wR * Fw
                # Moment -> moment
                f[2*idx+1] += wL * Mtheta
                f[2*(idx+1)+1] += wR * Mtheta
        
        for x_p, Fw in shear_pairs:
            if Fw != 0: add_load(float(x_p), float(Fw), 0.0)
        for x_p, Mth in moment_pairs:
            if Mth != 0: add_load(float(x_p), 0.0, float(Mth))
            
        return f

    # 3. BC Check
    def is_fixed(s: Support, k: int) -> bool:
        # k=0 -> w (transverse), k=1 -> theta (rotation)
        idx = trans_idx if k == 0 else rot_idx
        return s.type[idx] == "1"

    # --- Solve ---
    d, r = _solve_fem_1d(
        supports_sorted,
        dof_per_node=2,
        k_local_fn=k_hermite,
        f_global_fn=build_load_vector,
        is_fixed_fn=is_fixed
    )

    # Map reactions
    for i, s in enumerate(supports_sorted):
        w_dof = 2 * i
        t_dof = w_dof + 1
        s.reactions[shear_key] = float(r[w_dof])
        s.reactions[moment_key] = float(r[t_dof])

    return d



def get_all_loads(loads: LoadCase, beam: Beam1D) -> List[Tuple[float, str, float]]:
    """
    Returns a sorted list of all loads and support reactions as
        (x, type, magnitude)
    where repeated entries at the same x and of the same type are summed.
    """
    load_map: Dict[Tuple[float, str], float] = {}

    def add(x: float, type_: str, value: float):
        key = (float(x), type_)
        load_map[key] = load_map.get(key, 0.0) + float(value)

    # ---------------------------------------
    # Applied FORCES
    # ---------------------------------------
    for x, Fx in loads.Fxs:
        if Fx != 0:
            add(x, "Fx", Fx)

    for x, Fy in loads.Fys:
        if Fy != 0:
            add(x, "Fy", Fy)

    for x, Fz in loads.Fzs:
        if Fz != 0:
            add(x, "Fz", Fz)

    # ---------------------------------------
    # Applied MOMENTS
    # ---------------------------------------
    for x, Mx in loads.Mxs:
        if Mx != 0:
            add(x, "Mx", Mx)

    for x, My in loads.Mys:
        if My != 0:
            add(x, "My", My)

    for x, Mz in loads.Mzs:
        if Mz != 0:
            add(x, "Mz", Mz)

    # ---------------------------------------
    # Support reactions
    # ---------------------------------------
    for s in beam.supports:
        x = float(s.x)

        if s.reactions.get("Fx", 0.0) != 0.0: add(x, "Rx", s.reactions["Fx"])
        if s.reactions.get("Fy", 0.0) != 0.0: add(x, "Ry", s.reactions["Fy"])
        if s.reactions.get("Fz", 0.0) != 0.0: add(x, "Rz", s.reactions["Fz"])

        if s.reactions.get("Mx", 0.0) != 0.0: add(x, "RMx", s.reactions["Mx"])
        if s.reactions.get("My", 0.0) != 0.0: add(x, "RMy", s.reactions["My"])
        if s.reactions.get("Mz", 0.0) != 0.0: add(x, "RMz", s.reactions["Mz"])

    # ---------------------------------------
    # Convert to sorted list
    # ---------------------------------------
    out = [(x, t, val) for (x, t), val in load_map.items()]
    out.sort(key=lambda item: item[0])
    return out


# -------------------------------------------------
# Analysis Result Wrappers
# -------------------------------------------------
class Result:
    """
    Wraps analysis results (x, y) to provide convenient accessors.
    Behaves like a list of (x, y) tuples when iterated or printed.
    """
    def __init__(self, x: np.ndarray, values: np.ndarray):
        self._x = x
        self._values = values

    def __iter__(self):
        return zip(self._x.tolist(), self._values.tolist())

    def __repr__(self):
        return repr(list(self))

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return list(self)[idx]
        return (self._x[idx], self._values[idx])

    @property
    def max(self) -> float:
        return float(np.max(self._values))

    @property
    def min(self) -> float:
        return float(np.min(self._values))

    @property
    def mean(self) -> float:
        return float(np.mean(self._values))

    @property
    def range(self) -> float:
        return float(np.ptp(self._values))

    def at(self, x_loc: float) -> float:
        return float(np.interp(x_loc, self._x, self._values))


@dataclass
class AnalysisResult:
    _action: Result
    _stress: Result
    _displacement: Result

    @property
    def action(self) -> Result:
        """Force or Moment distribution (V, M, N, T)."""
        return self._action

    @property
    def stress(self) -> Result:
        """Stress distribution (sigma, tau)."""
        return self._stress

    @property
    def displacement(self) -> Result:
        """Displacement distribution (w, u, theta)."""
        return self._displacement


# -------------------------------------------------
# Internal Analysis Helpers
# -------------------------------------------------

def _accumulate_loads(
    points: np.ndarray,
    loads: List[Tuple[float, float]],
    moment_loads: Optional[List[Tuple[float, float]]] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Accumulate shear/axial forces and moments along the beam.
    
    Args:
        points: Array of x-coordinates to evaluate at.
        loads: List of (x, magnitude) for direct forces (Shear/Axial).
        moment_loads: List of (x, magnitude) for moment loads (Bending/Torsion).
                      If None, assumes no external moments contribute to the force summation directly.
                      
    Returns:
        Tuple[force_distribution, moment_distribution]
    """
    n_points = len(points)
    force_vals = np.zeros(n_points)
    moment_vals = np.zeros(n_points)
    
    has_moments = moment_loads is not None
    
    # Pre-sort for efficiency (though input should be sorted)
    # We'll iterate simply for clarity as n_points is usually small (<1000)
    # and number of loads is small.
    
    for i, x in enumerate(points):
        f_sum = 0.0
        m_sum = 0.0
        
        # Sum direct forces (V or N)
        # And their contribution to Moment (M = F * dist) if applicable
        for fx, fval in loads:
            if fx <= x + 1e-9:
                f_sum += fval
                # Moment arm contribution (only for transverse analysis)
                if has_moments:
                    m_sum += fval * (x - fx)
        
        # Sum direct moments (if any)
        if has_moments:
            if moment_loads is not None:
                for mx, mval in moment_loads:
                    if mx <= x + 1e-9:
                        m_sum += mval
        
        force_vals[i] = f_sum
        moment_vals[i] = m_sum
        
    return force_vals, moment_vals


def _hermite_displacement(
    points: np.ndarray,
    x_nodes: np.ndarray,
    d_vec: np.ndarray
) -> np.ndarray:
    """
    Interpolate displacements using Hermite shape functions.
    
    Args:
        points: Evaluation points.
        x_nodes: X-coordinates of the beam nodes (supports).
        d_vec: Global displacement vector [w1, theta1, w2, theta2, ...].
        
    Returns:
        Array of displacements at 'points'.
    """
    disps = np.zeros_like(points)
    
    for i, x in enumerate(points):
        # Find the element (node interval) containing x
        idx = np.searchsorted(x_nodes, x) - 1
        if idx < 0: idx = 0
        if idx >= len(x_nodes) - 1: idx = len(x_nodes) - 2
        
        xL = x_nodes[idx]
        xR = x_nodes[idx+1]
        L_elem = xR - xL
        
        if L_elem > 0:
            # Normalized coordinate xi in [0, 1]
            xi = (x - xL) / L_elem
            
            # Extract nodal displacements for this element
            # d_vec structure: [w0, t0, w1, t1, w2, t2, ...]
            dof_base = 2 * idx
            w1 = d_vec[dof_base]
            t1 = d_vec[dof_base + 1]
            w2 = d_vec[dof_base + 2]
            t2 = d_vec[dof_base + 3]
            
            # Hermite shape functions
            xi2 = xi * xi
            xi3 = xi2 * xi
            
            N1 = 1.0 - 3.0*xi2 + 2.0*xi3        # w1 influence
            N2 = L_elem * (xi - 2.0*xi2 + xi3)  # theta1 influence
            N3 = 3.0*xi2 - 2.0*xi3              # w2 influence
            N4 = L_elem * (-xi2 + xi3)          # theta2 influence
            
            disps[i] = N1*w1 + N2*t1 + N3*w2 + N4*t2
            
    return disps


def _linear_interpolation(
    points: np.ndarray,
    x_nodes: np.ndarray,
    d_vec: np.ndarray
) -> np.ndarray:
    """
    Interpolate axial/torsional displacements using linear shape functions.
    """
    return np.interp(points, x_nodes, d_vec)


# -------------------------------------------------
# Loaded beam wrapper
# -------------------------------------------------
from .plotter import plot_beam_diagram

@dataclass
class LoadedBeam:
    beam: Beam1D
    loads: LoadCase

    # State
    all_loads: List[Tuple[float, str, float]] = field(init=False)
    _d_y: np.ndarray = field(init=False, default=None)  # Displacements for axis 'y' (Mz bending)
    _d_z: np.ndarray = field(init=False, default=None)  # Displacements for axis 'z' (My bending)
    _d_x: np.ndarray = field(init=False, default=None)  # Displacements for axis 'x' (Axial)
    _d_rx: np.ndarray = field(init=False, default=None) # Displacements for axis 'rx' (Torsion)

    def __post_init__(self):
        """
        Build the full load set on instantiation:
        - convert distributed loads to point forces
        - solve for reactions (Fx, Fy/Fz, My/Mz)
        - assemble everything into all_loads
        """
        # 1) Reactions in each direction
        self._d_x, self._d_rx = solve_x_reactions(self.beam.supports, self.loads)
        self._d_y = solve_transverse_reactions(self.beam, self.loads, axis="y")
        self._d_z = solve_transverse_reactions(self.beam, self.loads, axis="z")

        # 2) Store combined applied loads + reactions
        self.all_loads = get_all_loads(self.loads, self.beam)

    def shear(self, axis: str, points: int = 100) -> AnalysisResult:
        return self._transverse_analysis(axis, points, mode="shear")

    def bending(self, axis: str, points: int = 100) -> AnalysisResult:
        return self._transverse_analysis(axis, points, mode="bending")

    def axial(self, points: int = 100) -> AnalysisResult:
        return self._axial_analysis(points, mode="axial")

    def torsion(self, points: int = 100) -> AnalysisResult:
        return self._axial_analysis(points, mode="torsion")

    def deflection(self, axis: str, points: int = 100) -> Result:
        return self._transverse_analysis(axis, points, mode="bending").displacement

    def von_mises(self, points: int = 100) -> Result:
        """
        Calculate the Von Mises stress distribution along the beam.
        
        Uses a conservative superposition of maximum stress components:
        sigma_vm = sqrt(sigma_max^2 + 3 * tau_max^2)
        
        where:
        sigma_max = |sigma_axial| + |sigma_bending_y| + |sigma_bending_z|
        tau_max = |tau_shear_y| + |tau_shear_z| + |tau_torsion|
        """
        # 1. Normal Stresses
        # Axial (Fx)
        r_ax = self.axial(points).stress
        
        # Bending about Z (Loads in Y) -> produces sigma_x varying with y
        r_by = self.bending("y", points).stress
        
        # Bending about Y (Loads in Z) -> produces sigma_x varying with z
        r_bz = self.bending("z", points).stress
        
        # Sum absolute maximums
        sigma_total = np.abs(r_ax._values) + np.abs(r_by._values) + np.abs(r_bz._values)
        
        # 2. Shear Stresses
        # Shear in Y
        r_sy = self.shear("y", points).stress
        
        # Shear in Z
        r_sz = self.shear("z", points).stress
        
        # Torsion
        r_tor = self.torsion(points).stress
        
        # Sum absolute maximums
        tau_total = np.abs(r_sy._values) + np.abs(r_sz._values) + np.abs(r_tor._values)
        
        # 3. Von Mises
        vm_values = np.sqrt(sigma_total**2 + 3 * tau_total**2)
        
        return Result(r_ax._x, vm_values)

    def plot(self, plot_stress: bool = False, plot_section: bool = True):
        """
        Plot the 3D beam diagram with loads.
        
        Args:
            plot_stress: If True, color the beam axis by von Mises stress
            plot_section: If True, draw the section outline at x=0
        """
        plot_beam_diagram(self, plot_stress=plot_stress, plot_section=plot_section)


    # ---------------------------------------------------------
    # Internal Analysis Logic
    # ---------------------------------------------------------
    def _transverse_analysis(self, axis: str, points: int, mode: str) -> AnalysisResult:
        """
        Calculate shear, bending, and displacement for transverse loading.
        """
        if axis not in ("y", "z"):
            raise ValueError("axis must be 'y' or 'z'")
            
        xs = np.linspace(0, self.beam.L, points)
        
        # Setup parameters based on axis
        if axis == "y":
            # Bending in x-y plane: Shear Fy, Moment Mz
            # Displacements: v (Uy), theta_z (Rz)
            shear_key, moment_key = "Fy", "Mz"
            d_vec = self._d_y
            I = self.beam.section.Iz
            c = self.beam.section.y_max
            A = self.beam.section.A
        else:
            # Bending in x-z plane: Shear Fz, Moment My
            # Displacements: w (Uz), theta_y (Ry)
            shear_key, moment_key = "Fz", "My"
            d_vec = self._d_z
            I = self.beam.section.Iy
            c = self.beam.section.z_max
            A = self.beam.section.A

        # 1. Filter Loads
        # Forces: Applied loads (F) + Reactions (R)
        # Moments: Applied moments (M) + Reactions (RM)
        relevant_forces = [
            (x, v) for x, t, v in self.all_loads 
            if t in (shear_key, f"R{shear_key[-1]}")
        ]
        relevant_moments = [
            (x, v) for x, t, v in self.all_loads 
            if t in (moment_key, f"RM{moment_key[-1]}")
        ]
        
        # 2. Calculate Internal Actions (V and M)
        V_vals, M_vals = _accumulate_loads(xs, relevant_forces, relevant_moments)

        # 3. Calculate Displacements
        # Reconstruct support x-locations for element interpolation
        supports_sorted = sorted(self.beam.supports, key=lambda s: s.x)
        x_nodes = np.array([s.x for s in supports_sorted])
        
        disps = np.zeros(points)
        if d_vec is not None:
            disps = _hermite_displacement(xs, x_nodes, d_vec)

        # 4. Package Results
        if mode == "shear":
            # Shear Force V and Shear Stress tau
            action_res = Result(xs, V_vals)
            stress_res = Result(xs, V_vals / A if A > 0 else np.zeros_like(V_vals))
            disp_res = Result(xs, disps)
            
        elif mode == "bending":
            # Bending Moment M and Bending Stress sigma
            action_res = Result(xs, M_vals)
            stress_val = (M_vals * c / I) if I > 0 else np.zeros_like(M_vals)
            stress_res = Result(xs, stress_val)
            disp_res = Result(xs, disps)
        else:
            raise ValueError(f"Unknown mode: {mode}")
            
        return AnalysisResult(action_res, stress_res, disp_res)

    def _axial_analysis(self, points: int, mode: str) -> AnalysisResult:
        """
        Calculate axial or torsional actions and displacements.
        """
        xs = np.linspace(0, self.beam.L, points)
        
        if mode == "axial":
            key = "Fx"
            prop_A = self.beam.section.A
            # stiffness = self.beam.material.E * prop_A # Not used in this path
            # dof_idx = 0 # Ux # Not used in this path
            # Radius/distance for stress calc (sigma = N/A) -> factor is 1/A
            stress_factor = 1.0 / prop_A if prop_A > 0 else 0.0
        else: # torsion
            key = "Mx"
            prop_J = self.beam.section.J
            # stiffness = self.beam.material.G * prop_J # Not used in this path
            # dof_idx = 3 # Rx # Not used in this path
            # Torsion stress tau = T*r/J. Use max dimension as conservative r
            r_max = max(abs(self.beam.section.y_max), abs(self.beam.section.z_max))
            stress_factor = r_max / prop_J if prop_J > 0 else 0.0
            
        # 1. Filter Loads
        reaction_key = "Rx" if key == "Fx" else "RMx"
        
        # For torsion, we specifically look for Mx and its reaction RMx
        # For axial, Fx and Rx
        relevant = [
            (x, v) for x, t, v in self.all_loads 
            if t in (key, reaction_key)
        ]
        
        # 2. Calculate Internal Actions (N or T)
        # Axial/Torsion doesn't have moment-arm accumulation from itself in this 1D view
        actions, _ = _accumulate_loads(xs, relevant, moment_loads=None)
        
        # 3. Calculate Displacements (u or theta)
        # Use the displacements computed by the FEM solver (linear interpolation)
        supports_sorted = sorted(self.beam.supports, key=lambda s: s.x)
        x_nodes = np.array([s.x for s in supports_sorted], dtype=float)
        
        if mode == "axial" and self._d_x is not None:
            disps = _linear_interpolation(xs, x_nodes, self._d_x)
        elif mode == "torsion" and self._d_rx is not None:
            disps = _linear_interpolation(xs, x_nodes, self._d_rx)
        else:
            # Fallback or empty
            disps = np.zeros_like(xs)
            
        # 4. Package Results
        action_res = Result(xs, actions)
        stress_res = Result(xs, actions * stress_factor)
        disp_res = Result(xs, disps)
        
        return AnalysisResult(action_res, stress_res, disp_res)

