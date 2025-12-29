from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

from beamy.analysis.assembly import AssemblyResult, assemble_global_system
from beamy.analysis.model import AnalysisModel
from beamy.analysis.recovery import element_axial_forces
from beamy.analysis.settings import AnalysisSettings
from beamy.analysis.solver_linear import LinearSolution, solve_linear_system
from beamy.loads.loadcase import LoadCase
from beamy.model.frame import Frame


@dataclass
class SecondOrderSolution:
    displacements: np.ndarray
    reactions: np.ndarray
    axial_forces: List[float]
    iterations: int
    converged: bool
    last_assembly: AssemblyResult


def solve_second_order(
    model: AnalysisModel,
    frame: Frame,
    loadcase: LoadCase,
    settings: AnalysisSettings,
    element_active: Optional[Sequence[bool]] = None,
) -> SecondOrderSolution:
    """
    Basic second-order (P-Delta / initial-stress) iteration.

    This is a fixed-point iteration on the geometric stiffness:
    - build Kg from axial forces from the previous iterate
    - solve (Km + Kg) u = F
    - repeat until u converges
    """
    total_dofs = len(model.nodes) * 6
    u_prev = np.zeros(total_dofs, dtype=float)
    converged = False
    last_assembly = assemble_global_system(model, frame, loadcase, element_active=element_active)

    for it in range(settings.max_iter):
        axial = element_axial_forces(model, u_prev)
        last_assembly = assemble_global_system(
            model,
            frame,
            loadcase,
            element_active=element_active,
            geometric_axial_forces=axial,
            include_geometric=True,
        )
        sol: LinearSolution = solve_linear_system(last_assembly)
        u_new = settings.relaxation * sol.displacements + (1.0 - settings.relaxation) * u_prev

        diff = np.linalg.norm(u_new - u_prev)
        denom = 1.0 + np.linalg.norm(u_new)
        if diff / denom < settings.tol:
            converged = True
            u_prev = u_new
            break
        u_prev = u_new

    axial_final = element_axial_forces(model, u_prev)
    reactions_full = last_assembly.k_full @ u_prev - last_assembly.f_full
    reactions = reactions_full[last_assembly.restrained_dofs] if last_assembly.restrained_dofs else np.array([])

    return SecondOrderSolution(
        displacements=u_prev,
        reactions=reactions,
        axial_forces=axial_final,
        iterations=it + 1,
        converged=converged,
        last_assembly=last_assembly,
    )

