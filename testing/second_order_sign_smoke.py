import numpy as np
from sectiony.library import rhs

from beamy import LoadCase, Material, FrameAnalysisSettings
from beamy.frame import FrameBuilder


def _find_node_id_at(frame, xyz, tol=1e-9):
    target = np.array(xyz, dtype=float)
    for nid, node in frame.nodes.items():
        if np.linalg.norm(node.position - target) <= tol:
            return nid
    raise RuntimeError(f"Could not find node at {xyz}")


def run() -> None:
    # Simple pinned-pinned column (along global X) with a mid-height lateral point load.
    # Compare second-order response under axial tension vs axial compression.
    L = 4.0

    steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
    sec = rhs(b=0.15, h=0.15, t=0.006, r=0.0)

    fb = FrameBuilder()
    fb.add("COL", (0.0, 0.0, 0.0), (L, 0.0, 0.0), sec, steel, orientation=(0, 1, 0))
    # Pinned-pinned for bending, but fix Rx at one end to remove the rigid-body torsional mode
    # (single isolated member otherwise has a free "spin" mode).
    fb.support_at((0.0, 0.0, 0.0), "111100")
    # Roller in axial (Ux) at the far end so axial load creates axial strain,
    # which is required for the current deformation-based axial force recovery.
    fb.support_at((L, 0.0, 0.0), "011000")
    frame = fb.build()

    n0 = _find_node_id_at(frame, (0.0, 0.0, 0.0))
    nL = _find_node_id_at(frame, (L, 0.0, 0.0))

    P = 200_000.0  # N
    H = 1_000.0  # N

    settings = FrameAnalysisSettings(
        analysis_method="SECOND_ORDER_ELASTIC",
        imperfection_model="none",
        max_iter=50,
        tol_u_rel=1e-8,
        tol_u_abs=1e-12,
        n_steps=5,
        relaxation_omega=0.8,
        cable_tension_only=False,
    )

    def solve_case(axial_right_fx: float) -> float:
        loads = LoadCase("case")
        loads.add_nodal_force(nL, np.array([axial_right_fx, 0.0, 0.0]), coords="global")
        loads.add_member_point_load(
            member_id="COL",
            position=0.5 * L,
            force=np.array([0.0, H, 0.0]),
            moment=np.array([0.0, 0.0, 0.0]),
            coords="global",
            position_type="absolute",
        )

        frame.analyze(loads, settings=settings)
        if not frame.analysis_result.converged:
            raise RuntimeError(f"Did not converge: {frame.analysis_result.warnings}")

        # Midpoint node should exist after load expansion.
        solve_state = frame._solve_state
        if solve_state is None:
            raise RuntimeError("Expected solve_state after analysis.")
        mid_id = _find_node_id_at(solve_state.expanded_frame, (0.5 * L, 0.0, 0.0), tol=1e-6)
        uy = float(frame.nodal_displacements[mid_id][1])
        return uy

    # Compression: apply axial force toward the fixed end at x=0.
    uy_comp = solve_case(-P)
    # Tension: apply axial force away from the fixed end.
    uy_tens = solve_case(+P)

    print(f"uy_comp = {uy_comp:.6e} m")
    print(f"uy_tens = {uy_tens:.6e} m")

    if abs(uy_comp) <= abs(uy_tens):
        raise AssertionError(
            "Expected compression case to be more flexible than tension case. "
            f"Got |uy_comp|={abs(uy_comp):.6e}, |uy_tens|={abs(uy_tens):.6e}"
        )


if __name__ == "__main__":
    run()
