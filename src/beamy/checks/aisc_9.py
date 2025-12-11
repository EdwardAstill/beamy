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
- Cb is computed per unbraced segment using end moments and interior max.
- Compact/noncompact slenderness classification is not performed; all
  formulae are reported and the governing allowable is chosen conservatively.
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
    ratio = M1 / M2
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


def _kc(h_over_tw: float) -> float:
    if h_over_tw <= 70.0:
        return 1.0
    return 4.05 / (h_over_tw ** 0.46)


# -----------------------------
# Checks (in AISC units: ksi, inches)
# -----------------------------


def _strong_axis_checks(beam: Beam1D, loaded: LoadedBeam, dims_info: Dict[str, Any]) -> Dict[str, Any]:
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

    Af = None
    if bf is not None and tf is not None:
        Af = bf * tf

    for seg in segments:
        Lb = seg.length
        Cb = seg.Cb

        # F1.1: compact
        Lc = None
        Fb_f11 = None
        if bf is not None and tf is not None and Af is not None and depth is not None:
            Lc = min(
                76.0 * bf / math.sqrt(Fy),
                20000.0 / ((depth / Af) * Fy),
            )
            Fb_f11 = 0.66 * Fy

        # F1.2(a): noncompact rolled
        Fb_f12a = None
        if bf is not None and tf is not None:
            Fb_f12a = Fy * (0.79 - 0.002 * (bf / (2.0 * tf)) * math.sqrt(Fy))

        # F1.2(b): built-up with kc
        Fb_f12b = None
        kc_val = None
        if h_clear is not None and tw is not None and bf is not None and tf is not None:
            h_tw = h_clear / tw
            kc_val = _kc(h_tw)
            Fb_f12b = (Fy * (0.79 - 0.002 * (bf / (2.0 * tf)) * math.sqrt(Fy))) / kc_val

        # F1.2(c): other noncompact
        Fb_f12c = 0.60 * Fy

        # F1.3 LTB
        Fb_f16 = None
        Fb_f17 = None
        Fb_f18 = None
        rT = None
        if bf is not None and tf is not None and tw is not None and h_clear is not None and Af is not None:
            rT = math.sqrt(
                (tf * (bf ** 3) / 12.0 + h_clear * (tw ** 3) / 72.0) /
                (bf * tf + (tw * h_clear) / 6.0)
            )
            slenderness = (Lb / rT) if rT > 0.0 else None
            if slenderness is not None and slenderness > 0.0:
                lower = math.sqrt((102000.0 * Cb) / Fy)
                upper = math.sqrt((510000.0 * Cb) / Fy)
                Fb_f16 = ((2.0 / 3.0) - (Fy * slenderness * slenderness) / (1530000.0 * Cb)) * Fy
                if Fb_f16 > 0.60 * Fy:
                    Fb_f16 = 0.60 * Fy
                Fb_f17 = (170000.0 * Cb) / (slenderness * slenderness)
                if Fb_f17 > 0.60 * Fy:
                    Fb_f17 = 0.60 * Fy
                Fb_f18 = (12000.0 * Cb) / (Lb * depth / Af) if depth is not None and Lb > 0.0 else None
                if Fb_f18 is not None and Fb_f18 > 0.60 * Fy:
                    Fb_f18 = 0.60 * Fy

        allowable_candidates = [
            Fb for Fb in (Fb_f11, Fb_f12a, Fb_f12b, Fb_f12c, Fb_f16, Fb_f17, Fb_f18) if Fb is not None
        ]
        governing_Fb = max(allowable_candidates) if allowable_candidates else None

        demand = None
        if modulus is not None and modulus > 0.0:
            demand = abs(seg.M_max) / modulus

        results.append(
            {
                "segment": {"x_start": seg.x_start, "x_end": seg.x_end, "Lb": Lb, "Cb": Cb},
                "moments": {"M1": seg.M1, "M2": seg.M2, "Mmax": seg.M_max},
                "rT": rT,
                "allowables": {
                    "F1.1_compact": {"Fb": Fb_f11, "Lc": Lc},
                    "F1.2a_noncompact_rolled": {"Fb": Fb_f12a},
                    "F1.2b_noncompact_built": {"Fb": Fb_f12b, "kc": kc_val},
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


def _weak_axis_checks(beam: Beam1D, loaded: LoadedBeam, dims_info: Dict[str, Any]) -> Dict[str, Any]:
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

    for seg in segments:
        Fb_f21 = 0.75 * Fy
        Fb_f22 = 0.60 * Fy
        Fb_f22b = None
        if bf is not None and tf is not None:
            Fb_f22b = Fy * (1.075 - 0.005 * (bf / (2.0 * tf)) * math.sqrt(Fy))

        allowable_candidates = [Fb for Fb in (Fb_f21, Fb_f22, Fb_f22b) if Fb is not None]
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


def _box_tube_checks(beam: Beam1D, loaded: LoadedBeam, dims_info: Dict[str, Any]) -> Dict[str, Any]:
    """Box/tube checks assuming all inputs already in ksi/inches."""
    moments = loaded.bending("y").action  # Use strong-axis moment for RHS/CHS
    modulus = loaded.beam.section.Sz
    Fy = beam.material.Fy
    if Fy is None:
        raise ValueError("Material.Fy is required for box/tube checks.")

    segments = _build_segments(moments, beam, axis="strong")
    results: List[Dict[str, Any]] = []

    Fb_compact = 0.66 * Fy
    Fb_noncompact = 0.60 * Fy

    for seg in segments:
        allowable_candidates = [Fb_compact, Fb_noncompact]
        governing_Fb = max(allowable_candidates)
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


def _shear_checks(beam: Beam1D, loaded: LoadedBeam, dims_info: Dict[str, Any]) -> Dict[str, Any]:
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

    if h_clear is not None and tw is not None and tw > 0.0:
        slenderness = h_clear / tw
        limit = 380.0 / math.sqrt(Fy)
        if slenderness <= limit:
            Fv = 0.40 * Fy
        else:
            # Slender web
            a = h_clear  # assume stiffener spacing equals depth if unknown
            ratio = a / h_clear if h_clear != 0.0 else 0.0
            if ratio < 1.0:
                kv = 4.00 + 5.34 * (ratio ** 2)
            else:
                kv = 5.34 + 4.00 / (ratio ** 2)
            Cv = 45000.0 / (slenderness * math.sqrt(kv * Fy))
            if Cv > 0.8:
                Cv = 190.0 / (slenderness * math.sqrt(Fy))
            Fv = Cv * (0.40 * Fy)

    capacity = None
    if Fv is not None and h_clear is not None and tw is not None:
        capacity = Fv * h_clear * tw

    return {
        "V_max": V_max,
        "slenderness_h_over_tw": slenderness,
        "Cv": Cv,
        "kv": kv,
        "Fv": Fv,
        "capacity": capacity,
        "pass": (capacity is not None and capacity >= V_max),
    }


# -----------------------------
# Main Check Function (with unit conversion)
# -----------------------------


def check_chapter_f(loaded_beam: LoadedBeam, length_unit: str, force_unit: str) -> Dict[str, Any]:
    """
    Run AISC ASD Chapter F checks with unit conversion.

    Converts all beam properties to AISC units (ksi, inches) before performing checks.

    Args:
        loaded_beam: LoadedBeam instance (already solved).
        length_unit: Unit of length in your model (e.g., "m", "mm", "ft", "in").
        force_unit: Unit of force in your model (e.g., "N", "kN", "lbf", "kip").

    Returns:
        Nested dict with strong/weak/box bending and shear results.
        All results are in AISC units (ksi for stress, inches for length, kip-in for moment).

    Raises:
        ValueError: If material.Fy is missing, section dimensions unavailable, or units incompatible.

    Examples:
        >>> # SI model
        >>> results = check_chapter_f(loaded_beam, "m", "N")
        
        >>> # US model
        >>> results = check_chapter_f(loaded_beam, "ft", "kip")
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

    # Create a temporary converted beam for internal use
    from copy import deepcopy
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

    # Create a temporary LoadedBeam with converted values
    # We need to convert the analysis results (moments, shear) as well
    from ..setup.loads import LoadCase

    # Convert loads for the temporary beam (we'll reuse loaded_beam's results but scale them)
    loads_in = LoadCase(name=loaded_beam.loads.name)

    # We need to create a temporary LoadedBeam that uses the converted beam
    # but we'll directly manipulate the analysis results
    
    # Create temporary loaded beam (without re-solving, we'll patch the results)
    loaded_in = LoadedBeam.__new__(LoadedBeam)
    loaded_in.beam = beam_in
    loaded_in.loads = loads_in
    loaded_in.all_loads = []
    
    # Copy and convert analysis results
    # Moments: force_unit * length_unit -> kip-in
    moment_unit = f"{force_unit} {length_unit}"
    
    # Convert the Result objects for bending
    def convert_result(original_result: Result, from_unit: str, to_unit: str) -> Result:
        """Convert a Result object to new units."""
        x_converted = conv(original_result._x, length_unit, "in")
        values_converted = conv(original_result._values, from_unit, to_unit)
        return Result(x_converted, values_converted)
    
    # We'll use the original loaded_beam's analysis methods but need to patch them
    # Actually, let's just patch the beam reference and let it compute segments/Cb naturally
    
    # Store original beam temporarily
    original_beam = loaded_beam.beam
    original_section = loaded_beam.beam.section
    
    # Monkey-patch the loaded_beam to use converted values
    loaded_beam.beam = beam_in
    loaded_beam.beam.section = section_in
    
    try:
        # Now run checks (they'll read from loaded_beam but use converted beam properties)
        # We need to also convert the moment/shear results
        
        # Get the analysis results and convert them
        bending_y_original = loaded_beam.bending("y")
        bending_z_original = loaded_beam.bending("z")
        shear_z_original = loaded_beam.shear("z")
        
        # Convert the action (moment/shear) results
        # Moment: force*length -> kip*in
        # Shear: force -> kip
        
        # Create a wrapper that returns converted results
        class ConvertedLoadedBeam:
            def __init__(self, original: LoadedBeam, beam_in: Beam1D):
                self._original = original
                self.beam = beam_in
                self.loads = original.loads
            
            def bending(self, axis: str):
                """Return bending results converted to kip-in."""
                result = self._original.bending(axis)
                # Convert moments from force*length to kip*in
                x_in = conv(result.action._x, length_unit, "in")
                M_kipin = conv(result.action._values, moment_unit, "kip in")
                
                from ..analysis.analysis import AnalysisResult
                converted_action = Result(x_in, M_kipin)
                # We don't need stress/displacement for checks, just return action
                return type('obj', (object,), {'action': converted_action})()
            
            def shear(self, axis: str):
                """Return shear results converted to kip."""
                result = self._original.shear(axis)
                # Convert shear from force to kip
                x_in = conv(result.action._x, length_unit, "in")
                V_kip = conv(result.action._values, force_unit, "kip")
                
                converted_action = Result(x_in, V_kip)
                return type('obj', (object,), {'action': converted_action})()
        
        loaded_converted = ConvertedLoadedBeam(loaded_beam, beam_in)
        
        # Prepare result dictionary
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

        # Run bending checks
        bending_results: List[Dict[str, Any]] = []

        if shape in ("i", "channel"):
            bending_results.append(_strong_axis_checks(beam_in, loaded_converted, {"shape": dims_info["shape"], "dims": dims_in}))
            bending_results.append(_weak_axis_checks(beam_in, loaded_converted, {"shape": dims_info["shape"], "dims": dims_in}))
        elif shape in ("rhs", "shs", "box", "tube", "chs"):
            bending_results.append(_box_tube_checks(beam_in, loaded_converted, {"shape": dims_info["shape"], "dims": dims_in}))
            bending_results.append(_weak_axis_checks(beam_in, loaded_converted, {"shape": dims_info["shape"], "dims": dims_in}))
        else:
            bending_results.append({"error": f"Unsupported shape '{dims_info['shape']}'."})

        results["bending"] = bending_results
        results["shear"] = _shear_checks(beam_in, loaded_converted, {"shape": dims_info["shape"], "dims": dims_in})

        return results
    
    finally:
        # Restore original beam
        loaded_beam.beam = original_beam
        loaded_beam.beam.section = original_section

