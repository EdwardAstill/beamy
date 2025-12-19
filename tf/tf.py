"""
TF frame model + plotting for the two loading conditions in tf/loading.md:

- Elevated:
  - Base RHS fixed
  - Apply force + moment to each post at applied height

- Lifted:
  - Two bottom nodes fixed in X/Y (Z free)
  - Lifting point constrained in Z (plus one extra DOF for numerical stability)
  - Apply force + moment to each post at applied height
  - Apply total base force (self-weight + sling) as distributed load on base RHS

All plots are saved as SVG into:
  - tf/outputs/elevated/
  - tf/outputs/lifted/
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
from sectiony.library import rhs, chs

from beamy import Material
from beamy.frame import Frame, FrameBuilder, FrameLoadCase, LoadedFrame


@dataclass(frozen=True)
class TfLoadSpec:
    name: str
    applied_height_m: float
    force_on_each_column_n: float
    moment_on_each_beam_nm: float
    base_force_n: float


MM_TO_M = 0.001


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _sum_lengths(frame: Frame, member_ids: Iterable[str]) -> float:
    total = 0.0
    for mid in member_ids:
        total += frame.member_lengths[mid]
    return total


def _print_reactions(loaded: LoadedFrame) -> None:
    total = np.zeros(3)
    for node_id, reaction in loaded.reactions.items():
        forces = reaction[:3]
        total += forces
        if float(np.linalg.norm(forces)) > 1e-6:
            print(f"Node {node_id}: FX={forces[0]: .2f} N, FY={forces[1]: .2f} N, FZ={forces[2]: .2f} N")
    print(f"Total reactions: FX={total[0]: .2f} N, FY={total[1]: .2f} N, FZ={total[2]: .2f} N")

def _frame_size_m(frame: Frame) -> float:
    pts = np.array(list(frame.node_positions.values()))
    spans = pts.max(axis=0) - pts.min(axis=0)
    return float(np.max(spans))


def _max_translation_displacement_m(loaded: LoadedFrame) -> float:
    return float(max(np.linalg.norm(d[:3]) for d in loaded.nodal_displacements.values()))


def _max_abs_component_displacement_m(loaded: LoadedFrame, component_idx: int) -> float:
    return float(max(abs(float(d[component_idx])) for d in loaded.nodal_displacements.values()))


def _auto_deflection_scale_factor(loaded: LoadedFrame, target_fraction_of_frame: float = 0.10) -> float:
    """
    Pick a plot scale factor so that max plotted deflection is ~10% of model size.

    This avoids the "everything is wildly warped" problem when the real displacements are a
    few mm but we use a large fixed scale factor.
    """
    frame_size = _frame_size_m(loaded.frame)
    max_u = _max_translation_displacement_m(loaded)
    if max_u < 1e-12:
        return 1.0
    target = target_fraction_of_frame * frame_size
    scale = target / max_u
    return float(np.clip(scale, 0.5, 200.0))


def _estimate_von_mises_range_pa(loaded: LoadedFrame, points_per_member: int = 15) -> tuple[float, float]:
    mn = float("inf")
    mx = 0.0
    for m in loaded.frame.members:
        res = loaded.get_member_results(m.id)
        for x in np.linspace(0.0, m.length, points_per_member):
            v = float(res.von_mises.at(float(x)))
            mn = min(mn, v)
            mx = max(mx, v)
    if not np.isfinite(mn):
        mn = 0.0
    if mx <= 0.0:
        mx = 1.0
    return mn, mx


def _mm(*coords: float) -> tuple[float, ...]:
    """Convert mm coordinates to meters."""
    return tuple(c * MM_TO_M for c in coords)


def build_frame_and_loads(
    spec: TfLoadSpec,
    support_mode: Literal["elevated", "lifted"],
    steel: Material,
    sling_mat: Material,
    base_section: object,
    post_section: object,
    sling_section: object,
    include_slings: bool,
) -> tuple[Frame, FrameLoadCase, dict[tuple[float, float, float], str]]:
    """
    Build a connected base + posts + slings frame and the matching loadcase.

    Uses FrameBuilder for cleaner coordinate-based member definition.

    Notes on supports:
    - Elevated: all Z=0 nodes are fully fixed (111111) to represent the fixed base RHS.
    - Lifted: two base nodes are fixed in X/Y only (110000), Z free. LIFT is fixed in Z, and
      we also fix LIFT in X to meet the minimum 6 constrained DOFs requirement.
    
    Returns:
        (Frame, FrameLoadCase, coord_to_node_id mapping)
    """
    z_load_mm = spec.applied_height_m / MM_TO_M  # convert back to mm for coordinate definitions
    z_300 = 300.0
    z_top = 1775.0
    if not (z_300 < z_load_mm < z_top):
        raise ValueError(f"Applied height must be between 300mm and 1775mm, got {z_load_mm:.1f} mm")

    fb = FrameBuilder()

    # === BASE FRAME (250x150x6 RHS) ===
    # BASE_L and BASE_R are split at intersections with BASE_C1 and BASE_C2
    # to create proper connectivity nodes at those crossing points
    fb.add("BASE_L_1", _mm(125, 0, 0), _mm(125, 750, 0), base_section, steel)
    fb.add("BASE_L_2", _mm(125, 750, 0), _mm(125, 1750, 0), base_section, steel)
    fb.add("BASE_L_3", _mm(125, 1750, 0), _mm(125, 2500, 0), base_section, steel)
    fb.add("BASE_R_1", _mm(1925, 0, 0), _mm(1925, 750, 0), base_section, steel)
    fb.add("BASE_R_2", _mm(1925, 750, 0), _mm(1925, 1750, 0), base_section, steel)
    fb.add("BASE_R_3", _mm(1925, 1750, 0), _mm(1925, 2500, 0), base_section, steel)
    fb.add("BASE_C1", _mm(125, 750, 0), _mm(1925, 750, 0), base_section, steel)
    fb.add("BASE_C2", _mm(125, 1750, 0), _mm(1925, 1750, 0), base_section, steel)

    # === POSTS (150x150x6 RHS) ===
    # Define each post as a single continuous member. The solver will auto-split
    # at z=300 (brace connection) and at the load height.
    post_ori = (1, 0, 0)  # local Y along global X
    post_locations = [(600, 750), (600, 1750), (1450, 750), (1450, 1750)]
    
    post_ids: list[str] = []
    for i, (px, py) in enumerate(post_locations, start=1):
        pid = f"P{i}"
        post_ids.append(pid)
        fb.add(pid, _mm(px, py, 0), _mm(px, py, 1775), post_section, steel, orientation=post_ori)

    # Connecting members between posts at Z=300 (with rotation releases)
    fb.add("BR1", _mm(600, 750, 300), _mm(600, 1750, 300), post_section, steel, releases="000111000111")
    fb.add("BR2", _mm(1450, 750, 300), _mm(1450, 1750, 300), post_section, steel, releases="000111000111")

    # === SLINGS (cables) ===
    if include_slings:
        lift_pt = _mm(1025, 1250, 2975)
        fb.add("SL1", _mm(600, 750, 1775), lift_pt, sling_section, sling_mat, element_type="cable")
        fb.add("SL2", _mm(600, 1750, 1775), lift_pt, sling_section, sling_mat, element_type="cable")
        fb.add("SL3", _mm(1450, 750, 1775), lift_pt, sling_section, sling_mat, element_type="cable")
        fb.add("SL4", _mm(1450, 1750, 1775), lift_pt, sling_section, sling_mat, element_type="cable")
        
        # Lift point support (lift point is a member end, so it is a real node)
        if support_mode == "lifted":
            fb.support_at(lift_pt, "101000")  # Ux + Uz fixed

    # === SUPPORTS ===
    if support_mode == "elevated":
        # Minimal supports to satisfy Frame validation; full base fixity is applied
        # via loads.support_member() after auto-splitting.
        fb.support_at(_mm(125, 0, 0), "111111")
        fb.support_at(_mm(1925, 2500, 0), "111111")
    else:  # lifted
        # Only two base nodes fixed in X/Y (Z free for lifting)
        fb.support_at(_mm(125, 0, 0), "110000")
        fb.support_at(_mm(1925, 2500, 0), "110000")

    # Build frame (coord mapping no longer needed)
    frame = fb.build()

    # === LOADS ===
    loads = FrameLoadCase(spec.name)
    
    # Apply forces and moments at load application points on each post (LOCAL coords)
    # Posts run along local x (global +Z), so a global -Z force becomes local Fx=-F.
    # The applied global Mx becomes a local My for these posts (given post_ori).
    post_rows = ["row_front", "row_back", "row_front", "row_back"]
    for pid, row in zip(post_ids, post_rows):
        moment_sign = 1.0 if row == "row_front" else -1.0
        loads.add_member_point_load(
            pid,
            position=spec.applied_height_m,
            force=np.array([-spec.force_on_each_column_n, 0.0, 0.0]),
            moment=np.array([0.0, moment_sign * spec.moment_on_each_beam_nm, 0.0]),
            coords="local",
            position_type="absolute",
        )

    # Elevated case: fix the entire base (all nodes created on those members)
    if support_mode == "elevated":
        for mid in ["BASE_L_1", "BASE_L_2", "BASE_L_3", "BASE_R_1", "BASE_R_2", "BASE_R_3", "BASE_C1", "BASE_C2"]:
            loads.support_member(mid, "111111")

    # Distributed self-weight on base members for lifted case
    if support_mode == "lifted":
        base_member_ids = ["BASE_L_1", "BASE_L_2", "BASE_L_3", "BASE_R_1", "BASE_R_2", "BASE_R_3", "BASE_C1", "BASE_C2"]
        base_total_length = _sum_lengths(frame, base_member_ids)
        w = spec.base_force_n / base_total_length  # N/m
        for mid in base_member_ids:
            m = frame.get_member(mid)
            w_global = np.array([0.0, 0.0, -w])
            w_local = m.transformation_matrix @ w_global
            loads.add_member_uniform_force(mid, w_local, coords="local")

    return frame, loads, {}


def run_case(
    output_dir: Path,
    spec: TfLoadSpec,
    support_mode: Literal["elevated", "lifted"],
    steel: Material,
    sling_mat: Material,
    base_section: object,
    post_section: object,
    sling_section: object,
    include_slings: bool,
) -> None:
    _ensure_dir(output_dir)
    frame, loads, _ = build_frame_and_loads(
        spec=spec,
        support_mode=support_mode,
        steel=steel,
        sling_mat=sling_mat,
        base_section=base_section,
        post_section=post_section,
        sling_section=sling_section,
        include_slings=include_slings,
    )

    print("=" * 70)
    print(f"{spec.name} ({support_mode})")
    print("=" * 70)
    print(f"Frame: {len(frame.nodes)} nodes, {len(frame.members)} members")
    print(loads)

    loaded = LoadedFrame(frame, loads)
    _print_reactions(loaded)

    # Extract P1 (first vertical post) as a loaded beam and analyze it
    print("\n" + "=" * 70)
    print("EXTRACTED BEAM ANALYSIS: Post P1")
    print("=" * 70)
    try:
        loaded_beam = loaded.to_loaded_beam("P1")
        print(f"Extracted P1 as LoadedBeam (L={loaded_beam.beam.L:.3f} m)")
        print(f"Loads: {loaded_beam.loads}")
        
        # Analyze deflection at midpoint and tip
        mid_x = loaded_beam.beam.L / 2.0
        tip_x = loaded_beam.beam.L
        
        # Get deflections from loaded beam analysis
        uy_result = loaded_beam.deflection("y")
        uz_result = loaded_beam.deflection("z")
        
        print(f"\nDeflections (1D beam analysis):")
        print(f"  Midspan (y): {uy_result.at(mid_x)*1000:.3f} mm")
        print(f"  Tip (y): {uy_result.at(tip_x)*1000:.3f} mm")
        print(f"  Midspan (z): {uz_result.at(mid_x)*1000:.3f} mm")
        print(f"  Tip (z): {uz_result.at(tip_x)*1000:.3f} mm")
        
        # Get bending moments
        by_result = loaded_beam.bending("y")
        bz_result = loaded_beam.bending("z")
        print(f"\nBending Moments (1D beam analysis):")
        print(f"  Midspan (My): {by_result.action.at(mid_x)/1e6:.3f} MNm")
        print(f"  Tip (My): {by_result.action.at(tip_x)/1e6:.3f} MNm")
        print(f"  Midspan (Mz): {bz_result.action.at(mid_x)/1e6:.3f} MNm")
        print(f"  Tip (Mz): {bz_result.action.at(tip_x)/1e6:.3f} MNm")
        
        # Get shear forces
        sy_result = loaded_beam.shear("y")
        sz_result = loaded_beam.shear("z")
        print(f"\nShear Forces (1D beam analysis):")
        print(f"  Start (Vy): {sy_result.action.at(1e-3)/1e3:.3f} kN")
        print(f"  Start (Vz): {sz_result.action.at(1e-3)/1e3:.3f} kN")
        
        # Compare with frame member analysis (use bundled segments if they exist)
        seg_ids = loaded._member_bundle.get("P1", ["P1"])
        if len(seg_ids) > 0:
            print(f"\nFrame Member (P1) segments: {seg_ids}")
            for seg_id in seg_ids:
                seg_res = loaded.get_member_results(seg_id)
                print(f"  {seg_id}: Max My = {seg_res.bending_y.abs_max/1e6:.3f} MNm, Max Mz = {seg_res.bending_z.abs_max/1e6:.3f} MNm")
        
        # Run AISC Chapter F check on the extracted beam
        print(f"\nRunning AISC Chapter F Check:")
        try:
            check_result = loaded_beam.check_aisc_chapter_f(length_unit="m", force_unit="N")
            # Print summary of check results
            print(f"  Overall Pass: {all(b.utilisation < 1.0 for b in check_result.bending) and check_result.shear.pass_}")
            for bending in check_result.bending:
                max_util = max((s.utilisation for s in bending.segments), default=0.0)
                print(f"  {bending.axis} axis utilization: {max_util:.3f}")
            print(f"  Shear utilization: {check_result.shear.utilisation:.3f} (V_max = {check_result.shear.V_max:.2f}, Capacity = {check_result.shear.capacity:.2f})")
        except Exception as check_err:
            print(f"  Check error: {check_err}")
        
    except Exception as e:
        print(f"Warning: Could not extract beam analysis: {e}")
        import traceback
        traceback.print_exc()

    # Make plotting robust/legible for both cases
    scale = _auto_deflection_scale_factor(loaded)
    vm_min, vm_max = _estimate_von_mises_range_pa(loaded)

    # Slings are element_type="cable", so tension-only behaviour is handled
    # inside LoadedFrame automatically.

    loaded.plot(deformed=True, scale_factor=scale, save_path=str(output_dir / "geometry.svg"))
    loaded.plot_deflection(
        scale_factor=scale,
        colormap="viridis",
        show_undeformed=True,
        save_path=str(output_dir / "deflection.svg"),
    )
    loaded.plot_von_mises(
        colormap="turbo",
        # Using steel.Fy here can make plots look "blank" when stresses are low.
        stress_limits=(0.0, 1.05 * vm_max),
        save_path=str(output_dir / "von_mises.svg"),
    )
    loaded.plot_aisc_utilization(save_path=str(output_dir / "aisc_utilization.png"))

    print(f"\nPlot scale: {scale:.2f}x (max disp = {_max_translation_displacement_m(loaded)*1000:.2f} mm)")
    print(f"Max |Ux| = {_max_abs_component_displacement_m(loaded, 0)*1000:.3f} mm")
    print(f"Von Mises range (Pa): min={vm_min:.3g}, max={vm_max:.3g}")
    print(f"Saved plots to: {output_dir}")


if __name__ == "__main__":
    # Materials + sections (SI units)
    steel = Material(name="Steel", E=200e9, G=80e9, Fy=345e6)
    sling_mat = Material(name="Sling", E=50e9, G=20e9, Fy=1000e6)

    base_section = rhs(b=0.150, h=0.250, t=0.006, r=0.0)
    post_section = rhs(b=0.150, h=0.150, t=0.006, r=0.0)
    sling_section = chs(d=0.020, t=0.002)

    # Values copied from tf/loading.md (units converted to SI where needed)
    elevated = TfLoadSpec(
        name="TF - Elevated",
        applied_height_m=1111.0 * MM_TO_M,
        force_on_each_column_n=6557.985000000001,
        moment_on_each_beam_nm=557428.7250000001 / 1000.0,  # Nmm -> Nm
        base_force_n=0.0,
    )
    lifted = TfLoadSpec(
        name="TF - Lifted",
        applied_height_m=875.0 * MM_TO_M,
        force_on_each_column_n=6557.985000000001,
        moment_on_each_beam_nm=557428.7250000001 / 1000.0,  # Nmm -> Nm
        base_force_n=16475.698800000002,
    )

    outputs_root = Path("tf") / "outputs"
    run_case(
        output_dir=outputs_root / "elevated",
        spec=elevated,
        support_mode="elevated",
        steel=steel,
        sling_mat=sling_mat,
        base_section=base_section,
        post_section=post_section,
        sling_section=sling_section,
        include_slings=False,  # elevated case doesn't use slings
    )
    run_case(
        output_dir=outputs_root / "lifted",
        spec=lifted,
        support_mode="lifted",
        steel=steel,
        sling_mat=sling_mat,
        base_section=base_section,
        post_section=post_section,
        sling_section=sling_section,
        include_slings=True,
    )
