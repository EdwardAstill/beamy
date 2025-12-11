"""
AISC 9th Edition Specification - Chapter F (ASD) checks.

This module evaluates a loaded beam (LoadedBeam) and returns a nested
dictionary with intermediate values and pass/fail summaries for:
    - F1 Strong-axis bending (I and channel shapes)
    - F2 Weak-axis bending
    - F3 Box / tube bending (RHS/CHS)
    - F4 Web shear

Notes / assumptions:
- Internally converts all values to AISC units (ksi, inches) before checks.
- Unbraced lengths come from support locations that restrain lateral
  translation and torsion:
    * Strong-axis (bending about z): braces require Uy and Rx fixed.
    * Weak-axis (bending about y): braces require Uz and Rx fixed.
- Cb is computed per unbraced segment using end moments and interior max, and
  is capped at 2.3 with Cb = 1.0 if an interior moment exceeds both ends.
- Basic compact/noncompact/slender classification is performed using B5.1-style
  limits for flanges/webs and box/tube/CHS walls. Slender elements are not
  treated with full AISC slender-element provisions; instead their allowable
  bending stress is limited to 0.60 Fy and a scope warning is issued.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from unity.core import conv
from ..analysis.analysis import LoadedBeam, Result
from ..setup.beam import Beam1D, Support

# -----------------------------
# Helpers
# -----------------------------


@dataclass
class BraceSegment:
    x_start: float
    x_end: float
    length: float
    Cb: float
    M1: float
    M2: float
    M_max: float


def _brace_positions(beam: Beam1D, axis: str) -> List[float]:
    """Find brace positions based on support DOFs."""
    brace_positions: List[float] = []
    for support in beam.supports:
        lateral_fixed = False
        torsion_fixed = support.type[3] == "1"  # Rx
        if axis == "strong":
            lateral_fixed = support.type[1] == "1"  # Uy
        else:
            lateral_fixed = support.type[2] == "1"  # Uz
        if lateral_fixed and torsion_fixed:
            brace_positions.append(float(support.x))
    brace_positions.append(0.0)
    brace_positions.append(float(beam.L))
    unique_positions = sorted(set(brace_positions))
    return unique_positions


def _cb_for_segment(moment_result: Result, x_start: float, x_end: float) -> Tuple[float, float, float, float]:
    """Compute Cb using segment end moments and interior maximum."""
    M1 = float(moment_result.at(x_start))
    M2 = float(moment_result.at(x_end))
    xs = np.linspace(x_start, x_end, 50)
    Ms = np.array([moment_result.at(float(xi)) for xi in xs])
    M_abs = np.abs(Ms)
    interior_max = float(np.max(M_abs))
    end_max = max(abs(M1), abs(M2))
    if interior_max > end_max + 1e-9:
        return 1.0, M1, M2, interior_max
    if abs(M2) < 1e-9:
        return 1.0, M1, M2, interior_max
    # AISC sign convention: Ratio is POSITIVE for reverse curvature (opposite signs).
    # Analysis moments have opposite signs for reverse curvature.
    # Therefore, we must invert the sign of the raw moment ratio.
    ratio = -1.0 * (M1 / M2)
    Cb = 1.75 + 1.05 * ratio + 0.30 * ratio * ratio
    if Cb > 2.3:
        Cb = 2.3
    if Cb < 0.0:
        Cb = 0.0
    return float(Cb), M1, M2, interior_max


def _build_segments(moment_result: Result, beam: Beam1D, axis: str) -> List[BraceSegment]:
    positions = _brace_positions(beam, axis)
    segments: List[BraceSegment] = []
    for i in range(len(positions) - 1):
        x_start = positions[i]
        x_end = positions[i + 1]
        Cb, M1, M2, Mmax = _cb_for_segment(moment_result, x_start, x_end)
        segments.append(
            BraceSegment(
                x_start=x_start,
                x_end=x_end,
                length=x_end - x_start,
                Cb=Cb,
                M1=M1,
                M2=M2,
                M_max=Mmax,
            )
        )
    return segments


def _infer_shape(dim_dict: Dict[str, float], section_name: str) -> str:
    keys = set(dim_dict.keys())
    name_lower = section_name.lower()
    if "tf" in keys and "tw" in keys:
        return "i"
    if "b" in keys and "h" in keys and "t" not in keys:
        return "solid_rect"
    if "d" in keys and "t" not in keys:
        return "solid_round"
    if "d" in keys and "t" in keys and len(keys) <= 3:
        return "chs"
    if "b" in keys and "h" in keys and "t" in keys:
        if "u " in name_lower or "channel" in name_lower:
            return "channel"
        if "rhs" in name_lower or "shs" in name_lower:
            return "rhs"
        return "rhs"
    return "unknown"


def _section_dims(beam: Beam1D) -> Dict[str, Any]:
    if not hasattr(beam.section, "dimensions"):
        raise ValueError("Section dimensions missing on section; expected 'dimensions' attribute.")
    dim_dict = beam.section.dimensions
    shape = _infer_shape(dim_dict, beam.section.name if hasattr(beam.section, "name") else "")
    return {"shape": shape, "dims": dim_dict}


def _dim(dim_dict: Dict[str, float], key: str) -> Optional[float]:
    if key in dim_dict:
        return dim_dict[key]
    return None


def _web_clear_depth(dim_dict: Dict[str, float]) -> Optional[float]:
    d_val = _dim(dim_dict, "d")
    h_val = _dim(dim_dict, "h")
    depth = d_val if d_val is not None else h_val
    tf_val = _dim(dim_dict, "tf")
    if depth is None:
        return None
    if tf_val is None:
        return depth
    return depth - 2.0 * tf_val


def _check_compactness(
    shape: str,
    dims: Dict[str, float],
    Fy: float,
) -> str:
    """
    Classify section as 'compact', 'noncompact', or 'slender' per the
    Table B5.1 limits provided in the Chapter F documentation excerpt.
    """
    status = "compact"

    bf = _dim(dims, "b")
    tf = _dim(dims, "tf")
    if tf is None:
        tf = _dim(dims, "t")
    tw = _dim(dims, "tw")
    if tw is None:
        tw = _dim(dims, "t")
    d = _dim(dims, "d")
    if d is None:
        d = _dim(dims, "h")
    h_clear = _web_clear_depth(dims)

    if shape in ("i", "channel") and bf and tf and tw and d:
        limit_p_f = 65.0 / math.sqrt(Fy)
        limit_r_f = 95.0 / math.sqrt(Fy)

        flange_ratio = bf / (2.0 * tf)
        if shape == "channel":
            flange_ratio = bf / tf

        if flange_ratio > limit_r_f:
            status = "slender"
        elif flange_ratio > limit_p_f:
            status = "noncompact"

        web_ratio = d / tw
        limit_p_w = 640.0 / math.sqrt(Fy)
        if web_ratio > limit_p_w:
            status = "slender"

    elif shape in ("rhs", "shs", "box", "tube") and bf and tf and d and tw:
        limit_p_f = 190.0 / math.sqrt(Fy)
        limit_r_f = 238.0 / math.sqrt(Fy)

        flange_ratio = (bf - 2.0 * tw) / tf if bf > 2.0 * tw else bf / tf
        if flange_ratio > limit_r_f:
            status = "slender"
        elif flange_ratio > limit_p_f:
            status = "noncompact"

        web_ratio = (d - 2.0 * tf) / tw if d > 2.0 * tf else d / tw
        limit_p_w = 640.0 / math.sqrt(Fy)
        if web_ratio > limit_p_w:
            status = "slender"
    elif shape == "chs" and d and tw:
        # CHS compactness per B5.1: compact if D/t <= 3300/Fy
        dt = d / tw
        limit_compact = 3300.0 / Fy
        if dt <= limit_compact:
            status = "compact"
        else:
            status = "slender"
    elif shape in ("solid_rect", "solid_round"):
        status = "compact"

    return status


# -----------------------------
# Checks (in AISC units: ksi, inches)
# -----------------------------


def _strong_axis_checks(
    beam: Beam1D,
    loaded: LoadedBeam,
    dims_info: Dict[str, Any],
    channel_major_axis: bool,
    classification_override: Optional[str],
) -> Dict[str, Any]:
    """Strong-axis checks assuming all inputs already in ksi/inches."""
    moments = loaded.bending("y").action  # Mz from Fy loads (strong axis)
    modulus = loaded.beam.section.Sz
    Fy = beam.material.Fy
    if Fy is None:
        raise ValueError("Material.Fy is required for strong-axis checks.")

    segments = _build_segments(moments, beam, axis="strong")
    results: List[Dict[str, Any]] = []

    dim_map = dims_info["dims"]
    bf = _dim(dim_map, "b")
    tf = _dim(dim_map, "tf") if _dim(dim_map, "tf") is not None else _dim(dim_map, "t")
    tw_val = _dim(dim_map, "tw")
    tw = tw_val if tw_val is not None else _dim(dim_map, "t")
    depth = _dim(dim_map, "d")
    if depth is None:
        depth = _dim(dim_map, "h")
    h_clear = _web_clear_depth(dim_map)

    # Check Compactness (B5.1)
    compactness = classification_override if classification_override else _check_compactness(
        dims_info["shape"],
        dim_map,
        Fy,
    )

    Af = None
    if bf is not None and tf is not None:
        Af = bf * tf

    for seg in segments:
        Lb = seg.length
        Cb = seg.Cb

        # --- Calculate Local Buckling / Yield Allowables ---
        
        # F1.1: Compact (0.66Fy)
        # Limit Lc is required to use this
        Lc = None
        Fb_f11 = None
        if bf is not None and tf is not None and Af is not None and depth is not None:
            Lc = min(
                76.0 * bf / math.sqrt(Fy),
                20000.0 / ((depth / Af) * Fy),
            )
            if compactness == "compact":
                Fb_f11 = 0.66 * Fy

        # F1.2: Noncompact (F1.2a / F1.2b / F1.2c)
        Fb_f12a = None
        Fb_f12c = 0.60 * Fy
        
        if compactness != "slender":
            if bf is not None and tf is not None:
                Fb_f12a = Fy * (0.79 - 0.002 * (bf / (2.0 * tf)) * math.sqrt(Fy))

        # Determine Local Allowable (ignoring LTB length for a moment)
        F_local = 0.60 * Fy # Fallback
        if compactness == "compact":
            F_local = 0.66 * Fy
        elif compactness == "noncompact":
             if Fb_f12a is not None:
                 F_local = Fb_f12a
        
        # --- Calculate LTB Allowables (F1.3) ---
        Fb_f16 = None
        Fb_f17 = None
        Fb_f18 = None
        rT = None
        F_LTB = 0.60 * Fy # Default max if LTB doesn't govern or calc fails
        
        if bf is not None and tf is not None and tw is not None and h_clear is not None and Af is not None:
            rT = math.sqrt(
                (tf * (bf ** 3) / 12.0 + h_clear * (tw ** 3) / 72.0) /
                (bf * tf + (tw * h_clear) / 6.0)
            )
            slenderness = (Lb / rT) if rT > 0.0 else None
            
            if slenderness is not None and slenderness > 0.0:
                lower = math.sqrt((102000.0 * Cb) / Fy)
                upper = math.sqrt((510000.0 * Cb) / Fy)

                # F1.3(a) Inelastic
                if (
                    slenderness >= lower
                    and slenderness <= upper
                    and not (dims_info["shape"] == "channel" and channel_major_axis)
                ):
                    Fb_f16 = ((2.0 / 3.0) - (Fy * slenderness * slenderness) / (1530000.0 * Cb)) * Fy
                    if Fb_f16 > 0.60 * Fy:
                        Fb_f16 = 0.60 * Fy

                # F1.3(b) Elastic
                if (slenderness >= upper or Fb_f16 is None) and not (dims_info["shape"] == "channel" and channel_major_axis):
                    Fb_f17 = (170000.0 * Cb) / (slenderness * slenderness)
                    if Fb_f17 > 0.60 * Fy:
                        Fb_f17 = 0.60 * Fy

                # F1.3(c) Local Flange Check (Eq F1-8)
                Fb_f18 = (12000.0 * Cb) / (Lb * depth / Af) if depth is not None and Lb > 0.0 else None
                if Fb_f18 is not None and Fb_f18 > 0.60 * Fy:
                    Fb_f18 = 0.60 * Fy
            
            channel_requires_f18 = dims_info["shape"] == "channel" and channel_major_axis
            if channel_requires_f18:
                F_LTB = Fb_f18 if Fb_f18 is not None else 0.60 * Fy
            else:
                ltb_main = Fb_f16 if Fb_f16 is not None else Fb_f17
                if ltb_main is None:
                    ltb_main = 0.60 * Fy
                ltb_flange = Fb_f18 if Fb_f18 is not None else 0.0
                F_LTB = max(ltb_main, ltb_flange)
                if F_LTB > 0.60 * Fy:
                    F_LTB = 0.60 * Fy

        # --- Determine Governing Fb ---
        governing_Fb = 0.60 * Fy
        
        if compactness == "compact" and Lc is not None and Lb <= Lc:
            governing_Fb = 0.66 * Fy
        else:
            if Lc is not None and Lb <= Lc and F_local is not None:
                governing_Fb = F_local
            else:
                governing_Fb = F_LTB
                if F_local is not None and F_local < governing_Fb:
                    governing_Fb = F_local

        demand = None
        if modulus is not None and modulus > 0.0:
            demand = abs(seg.M_max) / modulus

        results.append(
            {
                "segment": {"x_start": seg.x_start, "x_end": seg.x_end, "Lb": Lb, "Cb": Cb},
                "moments": {"M1": seg.M1, "M2": seg.M2, "Mmax": seg.M_max},
                "rT": rT,
                "classification": compactness,
                "allowables": {
                    "F1.1_compact": {"Fb": Fb_f11, "Lc": Lc},
                    "F1.2a_noncompact_rolled": {"Fb": Fb_f12a},
                    "F1.2c_other": {"Fb": Fb_f12c},
                    "F1.3a_inelastic": {"Fb": Fb_f16},
                    "F1.3b_elastic": {"Fb": Fb_f17},
                    "F1.3c_local_flange": {"Fb": Fb_f18},
                },
                "demand": {"Fb_required": demand},
                "pass": (demand is not None and governing_Fb is not None and governing_Fb >= demand),
                "governing_Fb": governing_Fb,
            }
        )

    return {"axis": "strong", "shape": dims_info["shape"], "segments": results}


def _weak_axis_checks(
    beam: Beam1D,
    loaded: LoadedBeam,
    dims_info: Dict[str, Any],
    classification_override: Optional[str],
) -> Dict[str, Any]:
    """Weak-axis checks assuming all inputs already in ksi/inches."""
    moments = loaded.bending("z").action  # My from Fz loads (weak axis)
    modulus = loaded.beam.section.Sy
    Fy = beam.material.Fy
    if Fy is None:
        raise ValueError("Material.Fy is required for weak-axis checks.")

    segments = _build_segments(moments, beam, axis="weak")
    results: List[Dict[str, Any]] = []

    dim_map = dims_info["dims"]
    bf = _dim(dim_map, "b")
    tf_candidate = _dim(dim_map, "tf")
    if tf_candidate is None:
        tf_candidate = _dim(dim_map, "t")
    tf = tf_candidate

    classification = classification_override if classification_override else _check_compactness(
        dims_info["shape"],
        dim_map,
        Fy,
    )

    for seg in segments:
        Fb_f21 = 0.75 * Fy
        Fb_f22 = 0.60 * Fy
        Fb_f22b = None
        if bf is not None and tf is not None:
            Fb_f22b = Fy * (1.075 - 0.005 * (bf / (2.0 * tf)) * math.sqrt(Fy))

        governing_Fb = None
        if classification == "compact":
            governing_Fb = Fb_f21
        else:
            allowable_candidates = [Fb_f22]
            if dims_info["shape"] == "i" and Fb_f22b is not None:
                allowable_candidates.append(Fb_f22b)
            governing_Fb = max(allowable_candidates) if allowable_candidates else None

        demand = None
        if modulus is not None and modulus > 0.0:
            demand = abs(seg.M_max) / modulus

        results.append(
            {
                "segment": {"x_start": seg.x_start, "x_end": seg.x_end, "Lb": seg.length, "Cb": seg.Cb},
                "moments": {"M1": seg.M1, "M2": seg.M2, "Mmax": seg.M_max},
                "allowables": {
                    "F2.1_compact": {"Fb": Fb_f21},
                    "F2.2_general": {"Fb": Fb_f22},
                    "F2.2b_noncompact_flange_I": {"Fb": Fb_f22b},
                },
                "demand": {"Fb_required": demand},
                "pass": (demand is not None and governing_Fb is not None and governing_Fb >= demand),
                "governing_Fb": governing_Fb,
            }
        )
    return {"axis": "weak", "shape": dims_info["shape"], "segments": results}


def _box_tube_checks(
    beam: Beam1D,
    loaded: LoadedBeam,
    dims_info: Dict[str, Any],
    classification_override: Optional[str],
) -> Dict[str, Any]:
    """Box/tube checks assuming all inputs already in ksi/inches."""
    moments = loaded.bending("y").action  # Use strong-axis moment for RHS/CHS
    modulus = loaded.beam.section.Sz
    Fy = beam.material.Fy
    if Fy is None:
        raise ValueError("Material.Fy is required for box/tube checks.")

    segments = _build_segments(moments, beam, axis="strong")
    results: List[Dict[str, Any]] = []

    dim_map = dims_info["dims"]
    compactness = classification_override if classification_override else _check_compactness(dims_info["shape"], dim_map, Fy)

    Fb_compact = 0.66 * Fy if compactness == "compact" else None
    Fb_noncompact = 0.60 * Fy if compactness != "slender" else 0.60 * Fy

    for seg in segments:
        limit_ok = True
        Lc_val = None
        b_val = _dim(dim_map, "b")
        h_val = _dim(dim_map, "h")
        t_val = _dim(dim_map, "t")
        tf_val = _dim(dim_map, "tf") if _dim(dim_map, "tf") is not None else t_val
        tw_val = _dim(dim_map, "tw") if _dim(dim_map, "tw") is not None else t_val
        if b_val is not None and h_val is not None and t_val is not None:
            depth_to_width_ok = h_val <= 6.0 * b_val
            thickness_ok = True
            if tf_val is not None and tw_val is not None:
                thickness_ok = tf_val <= 2.0 * tw_val
            limit_ok = depth_to_width_ok and thickness_ok
            ratio = 0.0
            if abs(seg.M2) > 0.0:
                ratio = seg.M1 / seg.M2
            Lc_val = (1950.0 + 1200.0 * ratio) / math.sqrt(Fy)
            lower_bound = 1200.0 * (b_val / math.sqrt(Fy))
            if Lc_val < lower_bound:
                Lc_val = lower_bound
            if seg.length > Lc_val:
                limit_ok = False

        allowable_candidates: List[float] = []
        if Fb_compact is not None and limit_ok:
            allowable_candidates.append(Fb_compact)
        if Fb_noncompact is not None:
            allowable_candidates.append(Fb_noncompact)
        governing_Fb = max(allowable_candidates) if allowable_candidates else None
        demand = None
        if modulus is not None and modulus > 0.0:
            demand = abs(seg.M_max) / modulus
        results.append(
            {
                "segment": {"x_start": seg.x_start, "x_end": seg.x_end, "Lb": seg.length, "Cb": seg.Cb},
                "moments": {"M1": seg.M1, "M2": seg.M2, "Mmax": seg.M_max},
                "allowables": {
                    "F3.1_compact": {"Fb": Fb_compact},
                    "F3.2_noncompact": {"Fb": Fb_noncompact},
                },
                "demand": {"Fb_required": demand},
                "pass": (demand is not None and governing_Fb >= demand),
                "governing_Fb": governing_Fb,
            }
        )
    return {"axis": "strong_box", "shape": dims_info["shape"], "segments": results}


def _shear_checks(
    beam: Beam1D,
    loaded: LoadedBeam,
    dims_info: Dict[str, Any],
    stiffener_spacing_in: Optional[float],
) -> Dict[str, Any]:
    """Shear checks assuming all inputs already in ksi/inches."""
    shear_res = loaded.shear("z").action  # Shear from vertical loads (Fz)
    V_max = float(np.max(np.abs(shear_res._values)))
    Fy = beam.material.Fy
    if Fy is None:
        raise ValueError("Material.Fy is required for shear checks.")

    dim_map = dims_info["dims"]
    tw_val = _dim(dim_map, "tw")
    tw = tw_val if tw_val is not None else _dim(dim_map, "t")
    h_clear = _web_clear_depth(dim_map)
    Fv = None
    Cv = None
    kv = None
    slenderness = None
    a_over_h = None

    if h_clear is not None and tw is not None and tw > 0.0:
        slenderness = h_clear / tw
        limit = 380.0 / math.sqrt(Fy)
        if slenderness <= limit:
            Fv = 0.40 * Fy
        else:
            if stiffener_spacing_in is not None and h_clear > 0.0:
                a_over_h = stiffener_spacing_in / h_clear
            if a_over_h is not None:
                if a_over_h < 1.0:
                    kv = 4.00 + 5.34 * (a_over_h * a_over_h)
                else:
                    kv = 5.34 + 4.00 / (a_over_h * a_over_h)
            else:
                kv = 5.34

            if kv is not None:
                Cv = 45000.0 / (slenderness * math.sqrt(kv * Fy))
                if Cv > 0.8:
                    Cv = 190.0 / (slenderness * math.sqrt(Fy))
                Fv = Cv * (0.40 * Fy)

    capacity = None
    if Fv is not None and h_clear is not None and tw is not None:
        capacity = Fv * h_clear * tw

    stiffeners_required = None
    if slenderness is not None and slenderness > 260.0 and capacity is not None:
        stiffeners_required = capacity < V_max

    return {
        "V_max": V_max,
        "slenderness_h_over_tw": slenderness,
        "Cv": Cv,
        "kv": kv,
        "Fv": Fv,
        "capacity": capacity,
        "a_over_h": a_over_h,
        "stiffeners_required": stiffeners_required,
        "pass": (capacity is not None and capacity >= V_max),
    }


# -----------------------------
# Main Check Function (with unit conversion)
# -----------------------------


def aisc_9_check(
    loaded_beam: LoadedBeam,
    length_unit: str,
    force_unit: str,
    *,
    channel_major_axis: bool = True,
    stiffener_spacing: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Run AISC ASD Chapter F checks with unit conversion.

    Converts beam properties to AISC units (ksi, inches) before performing checks.
    Additional parameters select spec branches:

    Args:
        loaded_beam: Solved LoadedBeam instance.
        length_unit: Length unit of the model (e.g., "m", "mm", "ft", "in").
        force_unit: Force unit of the model (e.g., "N", "kN", "lbf", "kip").
        channel_major_axis: If True, channels with Lb > Lc use Eq. F1-8 only.
        stiffener_spacing: Clear stiffener spacing (same units as beam) for F4 kv. If None, kv defaults to 5.34 (unstiffened).

    Returns:
        Nested dict with strong/weak/box bending and shear results in ksi/in units.

    Raises:
        ValueError: If required properties or units are missing/incompatible.
    """
    beam = loaded_beam.beam

    # Validate Fy exists
    if beam.material.Fy is None:
        raise ValueError("Material.Fy is required for AISC checks.")

    # Get section dimensions
    dims_info = _section_dims(beam)
    shape = dims_info["shape"].lower()

    # Build stress unit from force and length
    stress_unit_input = f"{force_unit} {length_unit}-2"

    # Convert Fy to ksi
    Fy_ksi = conv(beam.material.Fy, stress_unit_input, "kip in-2")

    # Convert E and G to ksi
    E_ksi = conv(beam.material.E, stress_unit_input, "kip in-2")
    G_ksi = conv(beam.material.G, stress_unit_input, "kip in-2")

    # Convert beam length to inches
    L_in = conv(beam.L, length_unit, "in")

    # Convert section properties
    A_in2 = conv(beam.section.A, f"{length_unit}2", "in2")
    Iy_in4 = conv(beam.section.Iy, f"{length_unit}4", "in4")
    Iz_in4 = conv(beam.section.Iz, f"{length_unit}4", "in4")
    J_in4 = conv(beam.section.J, f"{length_unit}4", "in4")
    
    Sy_in3 = None
    if hasattr(beam.section, "Sy") and beam.section.Sy is not None:
        Sy_in3 = conv(beam.section.Sy, f"{length_unit}3", "in3")
    
    Sz_in3 = None
    if hasattr(beam.section, "Sz") and beam.section.Sz is not None:
        Sz_in3 = conv(beam.section.Sz, f"{length_unit}3", "in3")

    # Convert section dimensions
    dims_in: Dict[str, float] = {}
    for key, val in dims_info["dims"].items():
        if val is not None:
            dims_in[key] = conv(val, length_unit, "in")

    # Convert support positions
    supports_in: List[Support] = []
    for support in beam.supports:
        x_in = conv(support.x, length_unit, "in")
        supports_in.append(Support(x=x_in, type=support.type))

    stiffener_spacing_in = conv(stiffener_spacing, length_unit, "in") if stiffener_spacing is not None else None

    # Create a temporary converted beam for internal use
    from ..setup.beam import Material, Beam1D
    from sectiony import Section

    # Create converted material
    material_in = Material(
        name=beam.material.name,
        E=E_ksi,
        G=G_ksi,
        Fy=Fy_ksi,
        transparency=beam.material.transparency,
    )

    # Create converted section
    section_in = Section(
        name=beam.section.name,
        A=A_in2,
        Iy=Iy_in4,
        Iz=Iz_in4,
        J=J_in4,
        Sy=Sy_in3,
        Sz=Sz_in3,
        y_max=conv(beam.section.y_max, length_unit, "in") if beam.section.y_max else None,
        z_max=conv(beam.section.z_max, length_unit, "in") if beam.section.z_max else None,
    )
    # Attach dimensions
    section_in.dimensions = dims_in

    # Create converted beam
    beam_in = Beam1D(
        L=L_in,
        material=material_in,
        section=section_in,
        supports=supports_in,
    )

    # Prepare unit strings for converting analysis results
    moment_unit = f"{force_unit} {length_unit}"
    
    # Store original beam temporarily
    original_beam = loaded_beam.beam
    original_section = loaded_beam.beam.section
    
    # Monkey-patch the loaded_beam to use converted values
    loaded_beam.beam = beam_in
    loaded_beam.beam.section = section_in
    
    try:
        class ConvertedLoadedBeam:
            def __init__(self, original: LoadedBeam, beam_in: Beam1D):
                self._original = original
                self.beam = beam_in
                self.loads = original.loads

            def bending(self, axis: str):
                result = self._original.bending(axis)
                x_in = conv(result.action._x, length_unit, "in")
                M_kipin = conv(result.action._values, moment_unit, "kip in")
                converted_action = Result(x_in, M_kipin)
                return type("obj", (object,), {"action": converted_action})()

            def shear(self, axis: str):
                result = self._original.shear(axis)
                x_in = conv(result.action._x, length_unit, "in")
                V_kip = conv(result.action._values, force_unit, "kip")
                converted_action = Result(x_in, V_kip)
                return type("obj", (object,), {"action": converted_action})()

        loaded_converted = ConvertedLoadedBeam(loaded_beam, beam_in)

        classification_default = _check_compactness(shape, dims_in, Fy_ksi)

        scope_warnings: List[str] = []
        if classification_default == "slender":
            scope_warnings.append(
                "One or more elements are classified as slender; slender-element behavior is approximated "
                "using noncompact (0.60 Fy) allowable stresses. Detailed slender provisions are not implemented."
            )

        h_clear_in = _web_clear_depth(dims_in)
        tw_in = _dim(dims_in, "tw")
        if tw_in is None:
            tw_in = _dim(dims_in, "t")
        if h_clear_in is not None and tw_in is not None and tw_in > 0.0:
            plate_girder_limit = 970.0 / math.sqrt(Fy_ksi)
            if (h_clear_in / tw_in) > plate_girder_limit:
                scope_warnings.append(
                    f"Web slenderness h/tw exceeds 970/sqrt(Fy) (value={h_clear_in/tw_in:.2f}); Chapter G plate girder rules may apply."
                )
        if Fy_ksi is not None and Fy_ksi > 65.0:
            scope_warnings.append("Fy exceeds 65 ksi; Chapter F limits may not apply for this material.")

        results: Dict[str, Any] = {
            "inputs": {
                "shape": dims_info["shape"],
                "original_units": {
                    "length": length_unit,
                    "force": force_unit,
                },
                "aisc_units": {
                    "length": "in",
                    "force": "kip",
                    "stress": "ksi",
                },
                "material_Fy_ksi": Fy_ksi,
                "E_ksi": E_ksi,
                "G_ksi": G_ksi,
                "section_props_in": {
                    "A_in2": A_in2,
                    "Iy_in4": Iy_in4,
                    "Iz_in4": Iz_in4,
                    "Sy_in3": Sy_in3,
                    "Sz_in3": Sz_in3,
                    "J_in4": J_in4,
                },
                "section_dims_in": dims_in,
            }
        }
        if scope_warnings:
            results["inputs"]["scope_warnings"] = scope_warnings

        bending_results: List[Dict[str, Any]] = []

        if shape in ("i", "channel"):
            bending_results.append(
                _strong_axis_checks(
                    beam_in,
                    loaded_converted,
                    {"shape": dims_info["shape"], "dims": dims_in},
                    channel_major_axis=channel_major_axis,
                    classification_override=classification_default,
                )
            )
            bending_results.append(
                _weak_axis_checks(
                    beam_in,
                    loaded_converted,
                    {"shape": dims_info["shape"], "dims": dims_in},
                    classification_override=classification_default,
                )
            )
        elif shape in ("solid_rect", "solid_round"):
            bending_results.append(
                _weak_axis_checks(
                    beam_in,
                    loaded_converted,
                    {"shape": dims_info["shape"], "dims": dims_in},
                    classification_override=classification_default,
                )
            )
        elif shape in ("rhs", "shs", "box", "tube", "chs"):
            bending_results.append(
                _box_tube_checks(
                    beam_in,
                    loaded_converted,
                    {"shape": dims_info["shape"], "dims": dims_in},
                    classification_override=classification_default,
                )
            )
            bending_results.append(
                _weak_axis_checks(
                    beam_in,
                    loaded_converted,
                    {"shape": dims_info["shape"], "dims": dims_in},
                    classification_override=classification_default,
                )
            )
        else:
            bending_results.append({"error": f"Unsupported shape '{dims_info['shape']}'."})

        results["bending"] = bending_results
        results["shear"] = _shear_checks(
            beam_in,
            loaded_converted,
            {"shape": dims_info["shape"], "dims": dims_in},
            stiffener_spacing_in=stiffener_spacing_in,
        )

        return results

    finally:
        loaded_beam.beam = original_beam
        loaded_beam.beam.section = original_section


def check_chapter_f(
    loaded_beam: LoadedBeam,
    length_unit: str,
    force_unit: str,
) -> Dict[str, Any]:
    """
    Backward-compatible wrapper for aisc_9_check using default spec choices.
    """
    return aisc_9_check(
        loaded_beam,
        length_unit,
        force_unit,
        channel_major_axis=True,
        stiffener_spacing=None,
    )

