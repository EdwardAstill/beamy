from __future__ import annotations

from beamy.analysis.assembly import assemble_global_system
from beamy.analysis.mesh import build_analysis_model
from beamy.analysis.recovery import recover_results
from beamy.analysis.settings import AnalysisSettings
from beamy.analysis.solver_linear import solve_linear_system, solve_tension_only
from beamy.analysis.solver_second_order import solve_second_order
from beamy.loads.loadcase import LoadCase
from beamy.model.frame import Frame
from beamy.results.frame_result import FrameResult


def run_analysis(frame: Frame, loadcase: LoadCase, settings: AnalysisSettings) -> FrameResult:
    settings.validate()
    model = build_analysis_model(frame, loadcase)

    has_cable = any(element.kind == "cable" for element in model.elements)

    if settings.solver == "linear":
        if has_cable:
            solution, active, assembly = solve_tension_only(model, frame, loadcase, settings)
            return recover_results(model, assembly, solution.displacements, loadcase.name)
        assembly = assemble_global_system(model, frame, loadcase)
        solution = solve_linear_system(assembly)
        return recover_results(model, assembly, solution.displacements, loadcase.name)

    # second order
    if has_cable:
        # approximate: first find cable active set from linear tension-only, then run second-order with that fixed.
        _, active, _ = solve_tension_only(model, frame, loadcase, settings)
        sol2 = solve_second_order(model, frame, loadcase, settings, element_active=active)
        return recover_results(model, sol2.last_assembly, sol2.displacements, loadcase.name)

    sol2 = solve_second_order(model, frame, loadcase, settings)
    return recover_results(model, sol2.last_assembly, sol2.displacements, loadcase.name)

