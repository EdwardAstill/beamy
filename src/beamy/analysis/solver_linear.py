from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from beamy.analysis.assembly import AssemblyResult, assemble_global_system
from beamy.analysis.model import AnalysisModel
from beamy.analysis.recovery import element_axial_forces
from beamy.analysis.settings import AnalysisSettings
from beamy.loads.loadcase import LoadCase
from beamy.model.frame import Frame


@dataclass
class LinearSolution:
    displacements: np.ndarray  # full vector
    reactions: np.ndarray
    free_dofs: List[int]
    restrained_dofs: List[int]


def solve_linear_system(assembly: AssemblyResult) -> LinearSolution:
    k_full = assembly.k_full
    f_full = assembly.f_full
    free = assembly.free_dofs
    restrained = assembly.restrained_dofs

    if not free:
        raise ValueError("no free degrees of freedom to solve")

    k_ff = k_full[np.ix_(free, free)]
    f_f = f_full[free]

    # known displacements at restrained dofs (supports default 0, prescribed override)
    u_full = np.zeros_like(f_full)
    for dof_idx, value in assembly.prescribed.items():
        u_full[int(dof_idx)] = float(value)

    if restrained:
        k_fc = k_full[np.ix_(free, restrained)]
        u_c = u_full[restrained]
        rhs = f_f - k_fc @ u_c
    else:
        rhs = f_f

    try:
        u_f = np.linalg.solve(k_ff, rhs)
    except np.linalg.LinAlgError:
        u_f, _, _, _ = np.linalg.lstsq(k_ff, rhs, rcond=None)

    for idx, value in zip(free, u_f):
        u_full[idx] = value

    reactions = k_full @ u_full - f_full
    reactions = reactions[restrained] if restrained else np.array([])

    return LinearSolution(
        displacements=u_full,
        reactions=reactions,
        free_dofs=free,
        restrained_dofs=restrained,
    )


def solve_tension_only(
    model: AnalysisModel,
    frame: Frame,
    loadcase: LoadCase,
    settings: AnalysisSettings,
    initial_active: Optional[Sequence[bool]] = None,
) -> Tuple[LinearSolution, List[bool], AssemblyResult]:
    """
    Simple tension-only active-set iteration for cable elements.

    - Treat cables as truss elements when active.
    - If axial force becomes compressive (N<0), deactivate the cable and re-solve.
    """
    active = list(initial_active) if initial_active is not None else [True for _ in model.elements]
    if len(active) != len(model.elements):
        raise ValueError("initial_active length must match number of elements")

    last_assembly = assemble_global_system(model, frame, loadcase, element_active=active)
    last_solution = solve_linear_system(last_assembly)

    for _ in range(settings.max_iter):
        axial = element_axial_forces(model, last_solution.displacements)
        changed = False
        for element in model.elements:
            if element.kind != "cable":
                continue
            should_be_active = axial[element.id] >= 0.0
            if active[element.id] != should_be_active:
                active[element.id] = should_be_active
                changed = True
        if not changed:
            return last_solution, active, last_assembly

        last_assembly = assemble_global_system(model, frame, loadcase, element_active=active)
        last_solution = solve_linear_system(last_assembly)

    return last_solution, active, last_assembly

