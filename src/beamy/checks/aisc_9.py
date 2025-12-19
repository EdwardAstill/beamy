"""AISC Specification checks (single-file implementation).

This module currently contains:
- Chapter F (ASD): Flexure + shear checks
- Chapter E (ASD): Axial compression (flexural buckling) checks
- Chapter H (ASD): Combined axial force and bending (interaction)
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union
import numpy as np
from unity.core import conv

from ..beam1d.analysis import LoadedMember
from ..beam1d.beam import Beam1D
from ..core.support import Support
from ..core.results import Result
from ..frame.results import MemberActionProfile

# -----------------------------------------------------------------------------
# 1. Result Structures (Dataclasses)
# -----------------------------------------------------------------------------

@dataclass
class BendingSegmentResult:
    x_start: float
    x_end: float
    Lb: float
    Cb: float
    M1: float
    M2: float
    Mmax: float
    classification: str
    governing_Fb: float
    demand_Fb: float
    pass_: bool
    allowables: Dict[str, Any]
    rT: Optional[float] = None

    @property
    def utilisation(self) -> float:
        if not self.governing_Fb: return 0.0
        return (self.demand_Fb or 0.0) / self.governing_Fb

@dataclass
class BendingResult:
    axis: str
    shape: str
    segments: List[BendingSegmentResult]

    @property
    def utilisation(self) -> float:
        return max((s.utilisation for s in self.segments), default=0.0)

    @property
    def pass_(self) -> bool:
        return all(s.pass_ for s in self.segments)

@dataclass
class ShearResult:
    V_max: float
    capacity: Optional[float]
    pass_: bool
    slenderness_h_over_tw: Optional[float] = None
    Cv: Optional[float] = None
    kv: Optional[float] = None
    Fv: Optional[float] = None
    a_over_h: Optional[float] = None
    stiffeners_required: Optional[bool] = None

    @property
    def utilisation(self) -> float:
        if not self.capacity: return 0.0
        return self.V_max / self.capacity

@dataclass
class AISC9Result:
    bending: List[BendingResult]
    shear: ShearResult
    inputs: Dict[str, Any]

    @property
    def utilisation(self) -> float:
        utils = [b.utilisation for b in self.bending] + [self.shear.utilisation]
        return max(utils) if utils else 0.0

    @property
    def pass_(self) -> bool:
        return all(b.pass_ for b in self.bending) and self.shear.pass_

    def info(self) -> Dict[str, Any]:
        """Returns a dictionary representation of the results for backward compatibility."""
        return {
            "inputs": self.inputs,
            "bending": [
                {
                    "axis": b.axis, "shape": b.shape,
                    "segments": [
                        {
                            "segment": {"x_start": s.x_start, "x_end": s.x_end, "Lb": s.Lb, "Cb": s.Cb},
                            "moments": {"M1": s.M1, "M2": s.M2, "Mmax": s.Mmax},
                            "rT": s.rT, "classification": s.classification,
                            "allowables": s.allowables, "demand": {"Fb_required": s.demand_Fb},
                            "pass": s.pass_, "governing_Fb": s.governing_Fb,
                        } for s in b.segments
                    ],
                } for b in self.bending
            ],
            "shear": {
                "V_max": self.shear.V_max, "slenderness_h_over_tw": self.shear.slenderness_h_over_tw,
                "Cv": self.shear.Cv, "kv": self.shear.kv, "Fv": self.shear.Fv,
                "capacity": self.shear.capacity, "a_over_h": self.shear.a_over_h,
                "stiffeners_required": self.shear.stiffeners_required, "pass": self.shear.pass_,
            },
        }


@dataclass
class CompressionSegmentResult:
    x_start: float
    x_end: float
    K: float
    L: float
    r_y: float
    r_z: float
    slenderness_y: float
    slenderness_z: float
    slenderness_governing: float
    Cc: float
    Fa: float
    P_comp_max: float
    P_allow: float

    @property
    def utilisation(self) -> float:
        if not self.P_allow:
            return 0.0
        return float(self.P_comp_max) / float(self.P_allow)

    @property
    def pass_(self) -> bool:
        return self.utilisation <= 1.0 + 1e-12


@dataclass
class AISCChapterEResult:
    segments: List[CompressionSegmentResult]
    inputs: Dict[str, Any]

    @property
    def utilisation(self) -> float:
        return max((s.utilisation for s in self.segments), default=0.0)

    @property
    def pass_(self) -> bool:
        return all(s.pass_ for s in self.segments)

    def info(self) -> Dict[str, Any]:
        return {
            "inputs": self.inputs,
            "compression": [
                {
                    "segment": {"x_start": s.x_start, "x_end": s.x_end, "L": s.L, "K": s.K},
                    "r": {"ry": s.r_y, "rz": s.r_z},
                    "slenderness": {
                        "KLr_y": s.slenderness_y,
                        "KLr_z": s.slenderness_z,
                        "governing": s.slenderness_governing,
                        "Cc": s.Cc,
                    },
                    "allowables": {"Fa": s.Fa, "P_allow": s.P_allow},
                    "demand": {"P_comp_max": s.P_comp_max},
                    "utilisation": s.utilisation,
                    "pass": s.pass_,
                }
                for s in self.segments
            ],
        }


@dataclass
class InteractionSegmentResult:
    """Result for a single segment under combined axial + bending (Chapter H)."""
    x_start: float
    x_end: float
    axial_type: str  # "compression" or "tension"
    fa: float  # computed axial stress
    fbx: float  # computed bending stress about strong axis
    fby: float  # computed bending stress about weak axis
    Fa: float  # allowable axial stress
    Fbx: float  # allowable bending stress about strong axis
    Fby: float  # allowable bending stress about weak axis
    F_prime_ex: Optional[float]  # Euler stress for x-axis bending plane
    F_prime_ey: Optional[float]  # Euler stress for y-axis bending plane
    Cmx: float  # moment gradient coefficient for x-axis
    Cmy: float  # moment gradient coefficient for y-axis
    H1_1: Optional[float]  # Buckling-modified interaction ratio (compression only)
    H1_2: Optional[float]  # Linear interaction ratio (compression only)
    H1_3: Optional[float]  # Simplified interaction ratio (compression, fa/Fa <= 0.15)
    H2_1: Optional[float]  # Tension + bending interaction ratio
    pass_: bool

    @property
    def utilisation(self) -> float:
        """Return the maximum interaction ratio."""
        ratios = [r for r in [self.H1_1, self.H1_2, self.H1_3, self.H2_1] if r is not None]
        return max(ratios) if ratios else 0.0


@dataclass
class AISCChapterHResult:
    """Result container for Chapter H combined stress checks."""
    segments: List[InteractionSegmentResult]
    inputs: Dict[str, Any]

    @property
    def utilisation(self) -> float:
        return max((s.utilisation for s in self.segments), default=0.0)

    @property
    def pass_(self) -> bool:
        return all(s.pass_ for s in self.segments)

    def info(self) -> Dict[str, Any]:
        return {
            "inputs": self.inputs,
            "interaction": [
                {
                    "segment": {"x_start": s.x_start, "x_end": s.x_end},
                    "axial_type": s.axial_type,
                    "stresses": {"fa": s.fa, "fbx": s.fbx, "fby": s.fby},
                    "allowables": {"Fa": s.Fa, "Fbx": s.Fbx, "Fby": s.Fby},
                    "euler": {"F_prime_ex": s.F_prime_ex, "F_prime_ey": s.F_prime_ey},
                    "Cm": {"Cmx": s.Cmx, "Cmy": s.Cmy},
                    "interaction_ratios": {
                        "H1-1": s.H1_1,
                        "H1-2": s.H1_2,
                        "H1-3": s.H1_3,
                        "H2-1": s.H2_1,
                    },
                    "utilisation": s.utilisation,
                    "pass": s.pass_,
                }
                for s in self.segments
            ],
        }

# -----------------------------------------------------------------------------
# 2. Engineering Logic (AISC Formulas in ksi/in)
# -----------------------------------------------------------------------------

def _dim(dim_dict: Dict[str, float], key: str) -> Optional[float]:
    return dim_dict.get(key)

def _web_clear_depth(dim_dict: Dict[str, float]) -> Optional[float]:
    depth = _dim(dim_dict, "d") or _dim(dim_dict, "h")
    tf_val = _dim(dim_dict, "tf")
    return (depth - 2.0 * tf_val) if (depth and tf_val) else depth

def _brace_positions(beam: Beam1D, axis: str) -> List[float]:
    pos = [0.0, float(beam.L)]
    for s in beam.supports:
        lat = s.type[1] == "1" if axis == "strong" else s.type[2] == "1"
        if lat and s.type[3] == "1": pos.append(float(s.x))
    return sorted(set(pos))

def _cb_for_segment(res: Result, x1: float, x2: float) -> Tuple[float, float, float, float]:
    m1, m2 = float(res.at(x1)), float(res.at(x2))
    mmax = float(np.max(np.abs([res.at(xi) for xi in np.linspace(x1, x2, 50)])))
    if mmax > max(abs(m1), abs(m2)) + 1e-9 or abs(m2) < 1e-9:
        return 1.0, m1, m2, mmax
    ratio = -1.0 * (m1 / m2)
    cb = min(2.3, max(0.0, 1.75 + 1.05 * ratio + 0.30 * ratio**2))
    return float(cb), m1, m2, mmax

def _check_compactness(shape: str, dims: Dict[str, float], Fy: float) -> str:
    status = "compact"
    bf, tf = _dim(dims, "b"), _dim(dims, "tf") or _dim(dims, "t")
    tw, d = _dim(dims, "tw") or _dim(dims, "t"), _dim(dims, "d") or _dim(dims, "h")
    if shape in ("i", "channel") and bf and tf and tw and d:
        fr = (bf / (2.0 * tf)) if shape == "i" else (bf / tf)
        if fr > 95.0 / math.sqrt(Fy): status = "slender"
        elif fr > 65.0 / math.sqrt(Fy): status = "noncompact"
        if d / tw > 640.0 / math.sqrt(Fy): status = "slender"
    elif shape in ("rhs", "shs", "box", "tube") and bf and tf and d and tw:
        fr = (bf - 2.0 * tw) / tf if bf > 2.0 * tw else bf / tf
        if fr > 238.0 / math.sqrt(Fy): status = "slender"
        elif fr > 190.0 / math.sqrt(Fy): status = "noncompact"
        if (d - 2.0 * tf) / tw > 640.0 / math.sqrt(Fy): status = "slender"
    elif shape == "chs" and d and tw:
        if d / tw > 3300.0 / Fy: status = "slender"
    return status

def _strong_axis_checks(
    length: float,
    moments: Result,
    dims_info: Dict[str, Any],
    Fy: float,
    modulus: Optional[float],
    brace_positions: List[float],
    channel_major: bool,
    classification: Optional[str],
) -> BendingResult:
    pos = sorted(brace_positions)
    results = []
    dim_map = dims_info["dims"]
    bf, tf = _dim(dim_map, "b"), _dim(dim_map, "tf") or _dim(dim_map, "t")
    tw, d = _dim(dim_map, "tw") or _dim(dim_map, "t"), _dim(dim_map, "d") or _dim(dim_map, "h")
    h_clear, Af = _web_clear_depth(dim_map), None
    if bf and tf:
        Af = bf * tf
    comp = classification or _check_compactness(dims_info["shape"], dim_map, Fy)

    for i in range(len(pos) - 1):
        x1, x2 = pos[i], pos[i + 1]
        cb, m1, m2, mmax = _cb_for_segment(moments, x1, x2)
        lb = x2 - x1
        lc = min(76.0 * bf / math.sqrt(Fy), 20000.0 / ((d / Af) * Fy)) if (bf and Af and d) else None

        f_local = 0.66 * Fy if comp == "compact" else 0.60 * Fy
        if comp == "noncompact" and bf and tf:
            f_local = Fy * (0.79 - 0.002 * (bf / (2.0 * tf)) * math.sqrt(Fy))

        f_ltb, rt = 0.60 * Fy, None
        if bf and tf and tw and h_clear and Af:
            rt = math.sqrt((tf * bf**3 / 12.0 + h_clear * tw**3 / 72.0) / (bf * tf + tw * h_clear / 6.0))
            slen = lb / rt if rt > 0 else 0
            if slen > 0:
                l, u = math.sqrt(102000.0 * cb / Fy), math.sqrt(510000.0 * cb / Fy)
                f16 = ((2 / 3) - (Fy * slen**2 / (1530000.0 * cb))) * Fy if l <= slen <= u else None
                f17 = (170000.0 * cb / slen**2) if slen >= u else None
                f18 = (12000.0 * cb / (lb * d / Af)) if (d and Af and lb > 0) else None
                cands = [f for f in [f16, f17, f18] if f is not None]
                if cands:
                    f_ltb = min(0.60 * Fy, max(cands))

        gov = f_local if (lc and lb <= lc) else min(f_local, f_ltb)
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(
            BendingSegmentResult(x1, x2, lb, cb, m1, m2, mmax, comp, gov, demand, gov >= demand, {}, rT=rt)
        )
    return BendingResult("strong", dims_info["shape"], results)


def _weak_axis_checks(
    length: float,
    moments: Result,
    dims_info: Dict[str, Any],
    Fy: float,
    modulus: Optional[float],
    brace_positions: List[float],
    classification: Optional[str],
) -> BendingResult:
    comp = classification or _check_compactness(dims_info["shape"], dims_info["dims"], Fy)
    results = []
    pos = sorted(brace_positions)
    for i in range(len(pos) - 1):
        x1, x2 = pos[i], pos[i + 1]
        cb, m1, m2, mmax = _cb_for_segment(moments, x1, x2)
        fb = 0.75 * Fy if comp == "compact" else 0.60 * Fy
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(BendingSegmentResult(x1, x2, x2 - x1, cb, m1, m2, mmax, comp, fb, demand, fb >= demand, {}))
    return BendingResult("weak", dims_info["shape"], results)


def _box_tube_checks(
    length: float,
    moments: Result,
    dims_info: Dict[str, Any],
    Fy: float,
    modulus: Optional[float],
    brace_positions: List[float],
    classification: Optional[str],
) -> BendingResult:
    comp = classification or _check_compactness(dims_info["shape"], dims_info["dims"], Fy)
    results = []
    pos = sorted(brace_positions)
    for i in range(len(pos) - 1):
        x1, x2 = pos[i], pos[i + 1]
        cb, m1, m2, mmax = _cb_for_segment(moments, x1, x2)
        fb = 0.66 * Fy if comp == "compact" else 0.60 * Fy
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(BendingSegmentResult(x1, x2, x2 - x1, cb, m1, m2, mmax, comp, fb, demand, fb >= demand, {}))
    return BendingResult("strong_box", dims_info["shape"], results)

def _shear_checks(shear_res: Result, dims_info: Dict[str, Any], Fy: float, stiffener_spacing_in: Optional[float]) -> ShearResult:
    vmax = float(np.max(np.abs(shear_res._values)))
    tw, h = _dim(dims_info["dims"], "tw") or _dim(dims_info["dims"], "t"), _web_clear_depth(dims_info["dims"])
    fv = 0.40 * Fy
    if h and tw and h / tw > 380 / math.sqrt(Fy):
        cv = 190.0 / ((h / tw) * math.sqrt(Fy))
        fv = cv * 0.40 * Fy
    cap = fv * h * tw if (h and tw) else 0.0
    return ShearResult(vmax, cap, cap >= vmax, Fv=fv, slenderness_h_over_tw=(h / tw if (h and tw) else None))


def _profile_from_loaded_beam(loaded_beam: LoadedMember) -> MemberActionProfile:
    beam = loaded_beam.beam
    try:
        torsion_action = loaded_beam.torsion().action
    except Exception:
        ax = loaded_beam.axial(points=3).action
        torsion_action = Result(ax._x, np.zeros_like(ax._values))

    return MemberActionProfile(
        member_id=getattr(beam, "id", "beam"),
        length=beam.L,
        material=beam.material,
        section=beam.section,
        axial=loaded_beam.axial(points=801).action,
        shear_y=loaded_beam.shear("y").action,
        shear_z=loaded_beam.shear("z").action,
        torsion=torsion_action,
        bending_y=loaded_beam.bending("y").action,
        bending_z=loaded_beam.bending("z").action,
    )

# -----------------------------------------------------------------------------
# 3. Main Interface & Unit Conversion
# -----------------------------------------------------------------------------

def aisc_9_check(
    subject: Union[LoadedMember, MemberActionProfile],
    length_unit: str,
    force_unit: str,
    *,
    channel_major_axis: bool = True,
    stiffener_spacing: Optional[float] = None,
    bracing_positions: Optional[Iterable[float]] = None,
) -> AISC9Result:
    profile = subject if isinstance(subject, MemberActionProfile) else _profile_from_loaded_beam(subject)
    if profile.material.Fy is None:
        raise ValueError("Material.Fy is required.")

    stress_in = f"{force_unit} {length_unit}-2"
    fy_ksi = conv(profile.material.Fy, stress_in, "kip in-2")
    e_ksi = conv(profile.material.E, stress_in, "kip in-2")
    g_ksi = conv(profile.material.G, stress_in, "kip in-2")
    length_in = conv(profile.length, length_unit, "in")

    brace_in = None
    if bracing_positions is not None:
        brace_in = [conv(x, length_unit, "in") for x in bracing_positions]
    elif not isinstance(subject, MemberActionProfile):
        # Fall back to supports on the beam for legacy path
        brace_in = _brace_positions(subject.beam, axis="strong")
        brace_in = [conv(x, length_unit, "in") for x in brace_in]
    strong_brace = _positions_with_ends(brace_in, length_in)
    weak_brace = strong_brace  # Use same positions unless user provides separate lists (not yet supported)

    dims_in = {k: conv(v, length_unit, "in") for k, v in profile.section.dimensions.items() if v is not None}
    shape = _infer_shape(dims_in, profile.section.name)

    # Convert action results
    axial_res = _convert_result_units(profile.axial, length_unit, "in", force_unit, "kip")
    shear_y_res = _convert_result_units(profile.shear_y, length_unit, "in", force_unit, "kip")
    shear_z_res = _convert_result_units(profile.shear_z, length_unit, "in", force_unit, "kip")
    torsion_res = _convert_result_units(profile.torsion, length_unit, "in", f"{force_unit} {length_unit}", "kip in")
    bending_y_res = _convert_result_units(profile.bending_y, length_unit, "in", f"{force_unit} {length_unit}", "kip in")
    bending_z_res = _convert_result_units(profile.bending_z, length_unit, "in", f"{force_unit} {length_unit}", "kip in")

    dims_info = {"shape": shape, "dims": dims_in}
    comp = _check_compactness(shape, dims_in, fy_ksi)
    bend_res: List[BendingResult] = []

    Sz_in = conv(profile.section.Sz, f"{length_unit}3", "in3") if profile.section.Sz else None
    Sy_in = conv(profile.section.Sy, f"{length_unit}3", "in3") if profile.section.Sy else None

    if shape in ("i", "channel"):
        bend_res.append(
            _strong_axis_checks(length_in, bending_y_res, dims_info, fy_ksi, Sz_in, strong_brace, channel_major_axis, comp)
        )
        bend_res.append(_weak_axis_checks(length_in, bending_z_res, dims_info, fy_ksi, Sy_in, weak_brace, comp))
    elif shape in ("rhs", "shs", "box", "tube", "chs"):
        bend_res.append(_box_tube_checks(length_in, bending_y_res, dims_info, fy_ksi, Sz_in, strong_brace, comp))
        bend_res.append(_weak_axis_checks(length_in, bending_z_res, dims_info, fy_ksi, Sy_in, weak_brace, comp))
    elif shape in ("solid_rect", "solid_round"):
        bend_res.append(_weak_axis_checks(length_in, bending_z_res, dims_info, fy_ksi, Sy_in, weak_brace, comp))

    stiff_in = conv(stiffener_spacing, length_unit, "in") if stiffener_spacing else None
    shear_res = _shear_checks(shear_z_res, dims_info, fy_ksi, stiff_in)

    return AISC9Result(
        bend_res,
        shear_res,
        {
            "shape": shape,
            "original_units": {"length": length_unit, "force": force_unit},
            "fy_ksi": fy_ksi,
            "E_ksi": e_ksi,
            "G_ksi": g_ksi,
        },
    )

def _infer_shape(dims: Dict[str, float], name: str) -> str:
    k, n = set(dims.keys()), name.lower()
    if "tf" in k and "tw" in k: return "i"
    if "b" in k and "h" in k: return "channel" if ("u " in n or "channel" in n) else ("rhs" if "t" in k else "solid_rect")
    if "d" in k: return "chs" if "t" in k else "solid_round"
    return "unknown"


def _positions_with_ends(unbraced_positions: Optional[Iterable[float]], L: float) -> List[float]:
    if unbraced_positions is None:
        return [0.0, float(L)]
    pos = sorted(set(float(x) for x in unbraced_positions))
    if not pos:
        return [0.0, float(L)]
    if pos[0] > 1e-9:
        pos = [0.0] + pos
    if pos[-1] < float(L) - 1e-9:
        pos = pos + [float(L)]
    out: List[float] = []
    for x in pos:
        out.append(min(max(0.0, x), float(L)))
    return sorted(set(out))


def _convert_result_units(res: Result, length_unit_from: str, length_unit_to: str, value_unit_from: str, value_unit_to: str) -> Result:
    return Result(
        conv(res._x, length_unit_from, length_unit_to),
        conv(res._values, value_unit_from, value_unit_to),
    )


def _fa_asd(E: float, Fy: float, slenderness: float) -> tuple[float, float]:
    """Return (Fa, Cc) for AISC ASD Chapter E.

    E and Fy in consistent stress units (here: ksi). Slenderness is KL/r.
    """
    if Fy <= 0.0 or E <= 0.0:
        return 0.0, 0.0
    Cc = math.sqrt((2.0 * math.pi**2 * E) / Fy)
    if slenderness < Cc:
        num = (1.0 - (slenderness**2) / (2.0 * Cc**2)) * Fy
        den = (5.0 / 3.0) + (3.0 * slenderness) / (8.0 * Cc) - (slenderness**3) / (8.0 * Cc**3)
        Fa = num / den if den != 0.0 else 0.0
        return float(Fa), float(Cc)
    Fa = (12.0 * math.pi**2 * E) / (23.0 * slenderness**2) if slenderness > 0.0 else 0.0
    return float(Fa), float(Cc)


def _f_prime_e_asd(E: float, K: float, Lb: float, r: float) -> float:
    """Compute Euler stress for Chapter H interaction (F'e).
    
    Per AISC H1.3:
    F'e = (12 π² E) / (23 (K Lb / r)²)
    
    All inputs in consistent units (ksi, inches).
    """
    if r <= 0.0 or Lb <= 0.0 or E <= 0.0:
        return 0.0
    slenderness = (K * Lb) / r
    if slenderness <= 0.0:
        return 0.0
    F_prime_e = (12.0 * math.pi**2 * E) / (23.0 * slenderness**2)
    return float(F_prime_e)


def _compute_cm(M1: float, M2: float, frame_type: str, has_transverse_loading: bool) -> float:
    """Compute moment gradient coefficient Cm per AISC H1.4.
    
    Args:
        M1: Smaller end moment (signed)
        M2: Larger end moment (signed)
        frame_type: "sidesway" or "braced"
        has_transverse_loading: True if member has transverse loads between ends
    
    Returns:
        Cm coefficient
    """
    if frame_type == "sidesway":
        return 0.85
    
    # Braced frame
    if has_transverse_loading:
        # Conservative: assume ends restrained against rotation
        return 0.85
    
    # No transverse loading, compute based on end moment ratio
    if abs(M2) < 1e-9:
        return 1.0
    
    ratio = M1 / M2
    Cm = 0.6 - 0.4 * ratio
    return float(Cm)


def aisc_chapter_h_check(
    subject: Union[LoadedMember, MemberActionProfile],
    length_unit: str,
    force_unit: str,
    *,
    Ky: float = 1.0,
    Kz: float = 1.0,
    unbraced_positions_y: Optional[Iterable[float]] = None,
    unbraced_positions_z: Optional[Iterable[float]] = None,
    Cmx: Optional[float] = None,
    Cmy: Optional[float] = None,
    frame_type: str = "braced",
    has_transverse_loading: bool = True,
    compression_sign: str = "negative",
) -> AISCChapterHResult:
    """ASD Chapter H: Combined axial force and bending (interaction).
    
    Args:
        subject: LoadedMember or MemberActionProfile with axial and bending results
        length_unit: Unit for lengths (e.g., "m", "ft", "in")
        force_unit: Unit for forces (e.g., "N", "kip", "lbf")
        Ky: Effective length factor for buckling about y-axis (strong axis)
        Kz: Effective length factor for buckling about z-axis (weak axis)
        unbraced_positions_y: Unbraced positions for strong-axis bending
        unbraced_positions_z: Unbraced positions for weak-axis bending
        Cmx: Moment gradient coefficient for strong axis (if None, computed)
        Cmy: Moment gradient coefficient for weak axis (if None, computed)
        frame_type: "sidesway" or "braced"
        has_transverse_loading: True if member has transverse loads
        compression_sign: "negative", "positive", or "abs"
    
    Returns:
        AISCChapterHResult with interaction check results
    """
    # Get profile
    profile = subject if isinstance(subject, MemberActionProfile) else _profile_from_loaded_beam(subject)
    if profile.material.Fy is None:
        raise ValueError("Material.Fy is required.")
    
    # Unit conversions
    stress_in = f"{force_unit} {length_unit}-2"
    fy_ksi = conv(profile.material.Fy, stress_in, "kip in-2")
    e_ksi = conv(profile.material.E, stress_in, "kip in-2")
    length_in = conv(profile.length, length_unit, "in")
    
    # Section properties
    A_in2 = conv(profile.section.A, f"{length_unit}2", "in2")
    Iy_in4 = conv(profile.section.Iy, f"{length_unit}4", "in4")
    Iz_in4 = conv(profile.section.Iz, f"{length_unit}4", "in4")
    Sz_in3 = conv(profile.section.Sz, f"{length_unit}3", "in3") if profile.section.Sz else None
    Sy_in3 = conv(profile.section.Sy, f"{length_unit}3", "in3") if profile.section.Sy else None
    
    r_y = math.sqrt(Iy_in4 / A_in2) if A_in2 > 0 and Iy_in4 > 0 else 0.0
    r_z = math.sqrt(Iz_in4 / A_in2) if A_in2 > 0 and Iz_in4 > 0 else 0.0
    
    # Convert action results to ksi/kip/in units
    axial_res = _convert_result_units(profile.axial, length_unit, "in", force_unit, "kip")
    bending_y_res = _convert_result_units(profile.bending_y, length_unit, "in", f"{force_unit} {length_unit}", "kip in")
    bending_z_res = _convert_result_units(profile.bending_z, length_unit, "in", f"{force_unit} {length_unit}", "kip in")
    
    # Get unbraced positions
    pos_y_in = _positions_with_ends(
        None if unbraced_positions_y is None else [conv(x, length_unit, "in") for x in unbraced_positions_y],
        length_in,
    )
    pos_z_in = _positions_with_ends(
        None if unbraced_positions_z is None else [conv(x, length_unit, "in") for x in unbraced_positions_z],
        length_in,
    )
    
    # Use the finer segmentation (more bracing points)
    all_positions = sorted(set(pos_y_in + pos_z_in))
    
    # Get dimensions and run Chapter F checks to get Fb values
    dims_in = {k: conv(v, length_unit, "in") for k, v in profile.section.dimensions.items() if v is not None}
    shape = _infer_shape(dims_in, profile.section.name)
    dims_info = {"shape": shape, "dims": dims_in}
    comp = _check_compactness(shape, dims_in, fy_ksi)
    
    # Get allowable bending stresses from Chapter F logic
    # For simplicity, we'll use conservative values based on compactness
    # In production, you'd want to run full Chapter F checks per segment
    if shape in ("i", "channel"):
        Fbx = 0.66 * fy_ksi if comp == "compact" else 0.60 * fy_ksi  # Strong axis
        Fby = 0.75 * fy_ksi if comp == "compact" else 0.60 * fy_ksi  # Weak axis
    elif shape in ("rhs", "shs", "box", "tube", "chs"):
        Fbx = 0.66 * fy_ksi if comp == "compact" else 0.60 * fy_ksi
        Fby = 0.75 * fy_ksi if comp == "compact" else 0.60 * fy_ksi
    else:
        Fbx = 0.60 * fy_ksi
        Fby = 0.60 * fy_ksi
    
    # Process segments
    segments: List[InteractionSegmentResult] = []
    
    for i in range(len(all_positions) - 1):
        x1, x2 = all_positions[i], all_positions[i + 1]
        L_seg = x2 - x1
        
        # Get axial force and moments at segment ends and max
        P_start = axial_res.at(x1)
        P_end = axial_res.at(x2)
        P_vals = [axial_res.at(xi) for xi in np.linspace(x1, x2, 20)]
        
        My_start = bending_y_res.at(x1)
        My_end = bending_y_res.at(x2)
        My_vals = [bending_y_res.at(xi) for xi in np.linspace(x1, x2, 20)]
        
        Mz_start = bending_z_res.at(x1)
        Mz_end = bending_z_res.at(x2)
        Mz_vals = [bending_z_res.at(xi) for xi in np.linspace(x1, x2, 20)]
        
        # Determine if compression or tension
        P_array = np.array(P_vals)
        if compression_sign == "positive":
            is_compression = np.any(P_array > 1e-6)
            P_max = float(np.max(P_array)) if is_compression else 0.0
        elif compression_sign == "abs":
            is_compression = True
            P_max = float(np.max(np.abs(P_array)))
        else:  # negative
            is_compression = np.any(P_array < -1e-6)
            P_max = float(np.max(-P_array)) if is_compression else 0.0
        
        # Compute stresses
        fa = P_max / A_in2 if A_in2 > 0 else 0.0
        
        My_max = float(np.max(np.abs(My_vals)))
        Mz_max = float(np.max(np.abs(Mz_vals)))
        
        fbx = My_max / Sz_in3 if Sz_in3 and Sz_in3 > 0 else 0.0  # Strong axis bending
        fby = Mz_max / Sy_in3 if Sy_in3 and Sy_in3 > 0 else 0.0  # Weak axis bending
        
        if is_compression:
            # Chapter H1: Compression + Bending
            
            # Get Fa from Chapter E formula
            # Find governing slenderness for this segment
            Lb_y = L_seg  # Strong axis unbraced length
            Lb_z = L_seg  # Weak axis unbraced length
            
            slenderness_y = (Ky * Lb_y / r_y) if r_y > 0 else float("inf")
            slenderness_z = (Kz * Lb_z / r_z) if r_z > 0 else float("inf")
            slenderness_max = max(slenderness_y, slenderness_z)
            
            if math.isfinite(slenderness_max):
                Fa, _ = _fa_asd(e_ksi, fy_ksi, slenderness_max)
            else:
                Fa = 0.0
            
            # Compute F'ex and F'ey (Euler stress for interaction)
            F_prime_ex = _f_prime_e_asd(e_ksi, Ky, Lb_y, r_y)
            F_prime_ey = _f_prime_e_asd(e_ksi, Kz, Lb_z, r_z)
            
            # Compute Cm values if not provided
            if Cmx is None:
                Cmx_val = _compute_cm(My_start, My_end, frame_type, has_transverse_loading)
            else:
                Cmx_val = Cmx
            
            if Cmy is None:
                Cmy_val = _compute_cm(Mz_start, Mz_end, frame_type, has_transverse_loading)
            else:
                Cmy_val = Cmy
            
            # Equation H1-1: Buckling-modified interaction
            term1 = fa / Fa if Fa > 0 else 0.0
            
            # Check for stability (prevent division by zero or negative denominators)
            denom_x = 1.0 - (fa / F_prime_ex) if F_prime_ex > 0 else 0.0
            denom_y = 1.0 - (fa / F_prime_ey) if F_prime_ey > 0 else 0.0
            
            if denom_x > 0.01 and Fbx > 0:
                term2 = (Cmx_val * fbx) / (denom_x * Fbx)
            else:
                term2 = float("inf") if fbx > 0 else 0.0
            
            if denom_y > 0.01 and Fby > 0:
                term3 = (Cmy_val * fby) / (denom_y * Fby)
            else:
                term3 = float("inf") if fby > 0 else 0.0
            
            H1_1_val = term1 + term2 + term3 if math.isfinite(term2) and math.isfinite(term3) else float("inf")
            
            # Equation H1-2: Linear interaction
            H1_2_val = (fa / (0.60 * fy_ksi)) + (fbx / Fbx if Fbx > 0 else 0.0) + (fby / Fby if Fby > 0 else 0.0)
            
            # Equation H1-3: Simplified (if fa/Fa <= 0.15)
            H1_3_val = None
            if Fa > 0 and (fa / Fa) <= 0.15:
                H1_3_val = (fa / Fa) + (fbx / Fbx if Fbx > 0 else 0.0) + (fby / Fby if Fby > 0 else 0.0)
            
            # Check passes if all applicable equations <= 1.0
            checks = [H1_1_val, H1_2_val]
            if H1_3_val is not None:
                checks = [H1_3_val]  # Can use simplified equation instead
            
            pass_check = all(c <= 1.0 + 1e-9 for c in checks if math.isfinite(c))
            
            segments.append(InteractionSegmentResult(
                x_start=conv(x1, "in", length_unit),
                x_end=conv(x2, "in", length_unit),
                axial_type="compression",
                fa=fa,
                fbx=fbx,
                fby=fby,
                Fa=Fa,
                Fbx=Fbx,
                Fby=Fby,
                F_prime_ex=F_prime_ex,
                F_prime_ey=F_prime_ey,
                Cmx=Cmx_val,
                Cmy=Cmy_val,
                H1_1=H1_1_val if math.isfinite(H1_1_val) else None,
                H1_2=H1_2_val,
                H1_3=H1_3_val,
                H2_1=None,
                pass_=pass_check,
            ))
        
        else:
            # Chapter H2: Tension + Bending
            # Get allowable tensile stress (0.60 Fy for ASD)
            Ft = 0.60 * fy_ksi
            
            # Equation H2-1
            H2_1_val = (fa / Ft if Ft > 0 else 0.0) + (fbx / Fbx if Fbx > 0 else 0.0) + (fby / Fby if Fby > 0 else 0.0)
            
            pass_check = H2_1_val <= 1.0 + 1e-9
            
            segments.append(InteractionSegmentResult(
                x_start=conv(x1, "in", length_unit),
                x_end=conv(x2, "in", length_unit),
                axial_type="tension",
                fa=fa,
                fbx=fbx,
                fby=fby,
                Fa=Ft,  # Store tension allowable in Fa field
                Fbx=Fbx,
                Fby=Fby,
                F_prime_ex=None,
                F_prime_ey=None,
                Cmx=0.0,
                Cmy=0.0,
                H1_1=None,
                H1_2=None,
                H1_3=None,
                H2_1=H2_1_val,
                pass_=pass_check,
            ))
    
    return AISCChapterHResult(
        segments=segments,
        inputs={
            "Ky": Ky,
            "Kz": Kz,
            "Cmx": Cmx,
            "Cmy": Cmy,
            "frame_type": frame_type,
            "has_transverse_loading": has_transverse_loading,
            "compression_sign": compression_sign,
            "original_units": {"length": length_unit, "force": force_unit},
            "fy_ksi": fy_ksi,
            "E_ksi": e_ksi,
        },
    )


def aisc_chapter_e_check(
    loaded_beam: Union[LoadedMember, MemberActionProfile],
    length_unit: str,
    force_unit: str,
    *,
    K: float = 1.0,
    unbraced_positions: Optional[Iterable[float]] = None,
    compression_sign: str = "negative",
) -> AISCChapterEResult:
    """ASD Chapter E: axial compression (flexural buckling).

    This is axial-only. For combined P+M interaction, use Chapter H.
    """
    beam = loaded_beam.beam if isinstance(loaded_beam, LoadedMember) else None
    material = beam.material if beam else loaded_beam.material
    section = beam.section if beam else loaded_beam.section
    length_native = beam.L if beam else loaded_beam.length
    if material.Fy is None:
        raise ValueError("Material.Fy is required.")

    stress_in = f"{force_unit} {length_unit}-2"
    fy_ksi = conv(material.Fy, stress_in, "kip in-2")
    e_ksi = conv(material.E, stress_in, "kip in-2")

    A_in2 = conv(section.A, f"{length_unit}2", "in2")
    Iy_in4 = conv(section.Iy, f"{length_unit}4", "in4")
    Iz_in4 = conv(section.Iz, f"{length_unit}4", "in4")

    r_y = math.sqrt(Iy_in4 / A_in2) if A_in2 > 0 and Iy_in4 > 0 else 0.0
    r_z = math.sqrt(Iz_in4 / A_in2) if A_in2 > 0 and Iz_in4 > 0 else 0.0

    if isinstance(loaded_beam, MemberActionProfile):
        ax = _convert_result_units(loaded_beam.axial, length_unit, "in", force_unit, "kip")
    else:
        ax = loaded_beam.axial(points=801).action
    p_kip = np.asarray(conv(ax._values, force_unit, "kip"), dtype=float)
    if compression_sign == "positive":
        p_comp = np.maximum(p_kip, 0.0)
    elif compression_sign == "abs":
        p_comp = np.abs(p_kip)
    else:
        p_comp = np.maximum(-p_kip, 0.0)
    p_comp_max = float(np.max(p_comp)) if p_comp.size else 0.0

    L_in = conv(length_native, length_unit, "in")
    pos_in = _positions_with_ends(
        None if unbraced_positions is None else [conv(x, length_unit, "in") for x in unbraced_positions],
        L_in,
    )

    segs: List[CompressionSegmentResult] = []
    for x1, x2 in zip(pos_in[:-1], pos_in[1:]):
        L_seg = float(x2 - x1)
        slender_y = (float(K) * L_seg / r_y) if r_y > 0 else float("inf")
        slender_z = (float(K) * L_seg / r_z) if r_z > 0 else float("inf")
        slender_g = max(slender_y, slender_z)
        Fa, Cc = _fa_asd(float(e_ksi), float(fy_ksi), slender_g if math.isfinite(slender_g) else 0.0)
        P_allow = float(Fa) * float(A_in2) if A_in2 > 0 else 0.0
        segs.append(
            CompressionSegmentResult(
                x_start=float(x1),
                x_end=float(x2),
                K=float(K),
                L=L_seg,
                r_y=float(r_y),
                r_z=float(r_z),
                slenderness_y=float(slender_y),
                slenderness_z=float(slender_z),
                slenderness_governing=float(slender_g),
                Cc=float(Cc),
                Fa=float(Fa),
                P_comp_max=float(p_comp_max),
                P_allow=float(P_allow),
            )
        )

    return AISCChapterEResult(
        segments=segs,
        inputs={
            "K": float(K),
            "compression_sign": compression_sign,
            "original_units": {"length": length_unit, "force": force_unit},
            "fy_ksi": float(fy_ksi),
            "E_ksi": float(e_ksi),
        },
    )

def run(loaded_beam: LoadedMember, **kwargs) -> Any:
    """Dispatch to Chapter F (default), Chapter E, or Chapter H.

    - `chapter="f"` (default): returns `AISC9Result`
    - `chapter="e"`: returns `AISCChapterEResult`
    - `chapter="h"`: returns `AISCChapterHResult`
    """
    chapter = str(kwargs.get("chapter", "f")).strip().lower()
    if chapter in ("f", "chapter_f", "chapterf"):
        return aisc_9_check(
            loaded_beam,
            kwargs["length_unit"],
            kwargs["force_unit"],
            channel_major_axis=kwargs.get("channel_major_axis", True),
            stiffener_spacing=kwargs.get("stiffener_spacing"),
        )
    if chapter in ("e", "chapter_e", "chaptere"):
        return aisc_chapter_e_check(
            loaded_beam,
            kwargs["length_unit"],
            kwargs["force_unit"],
            K=kwargs.get("K", 1.0),
            unbraced_positions=kwargs.get("unbraced_positions"),
            compression_sign=kwargs.get("compression_sign", "negative"),
        )
    if chapter in ("h", "chapter_h", "chapterh"):
        return aisc_chapter_h_check(
            loaded_beam,
            kwargs["length_unit"],
            kwargs["force_unit"],
            Ky=kwargs.get("Ky", 1.0),
            Kz=kwargs.get("Kz", 1.0),
            unbraced_positions_y=kwargs.get("unbraced_positions_y"),
            unbraced_positions_z=kwargs.get("unbraced_positions_z"),
            Cmx=kwargs.get("Cmx"),
            Cmy=kwargs.get("Cmy"),
            frame_type=kwargs.get("frame_type", "braced"),
            has_transverse_loading=kwargs.get("has_transverse_loading", True),
            compression_sign=kwargs.get("compression_sign", "negative"),
        )
    raise ValueError(f"Unknown AISC chapter: {chapter!r}")

def check_chapter_f(loaded_beam: LoadedMember, l_u: str, f_u: str) -> Dict[str, Any]:
    return aisc_9_check(loaded_beam, l_u, f_u).info()


def check_chapter_e(
    loaded_beam: LoadedMember,
    l_u: str,
    f_u: str,
    *,
    K: float = 1.0,
    unbraced_positions: Optional[Iterable[float]] = None,
    compression_sign: str = "negative",
) -> Dict[str, Any]:
    return aisc_chapter_e_check(
        loaded_beam,
        l_u,
        f_u,
        K=K,
        unbraced_positions=unbraced_positions,
        compression_sign=compression_sign,
    ).info()


def check_chapter_h(
    loaded_beam: LoadedMember,
    l_u: str,
    f_u: str,
    *,
    Ky: float = 1.0,
    Kz: float = 1.0,
    unbraced_positions_y: Optional[Iterable[float]] = None,
    unbraced_positions_z: Optional[Iterable[float]] = None,
    Cmx: Optional[float] = None,
    Cmy: Optional[float] = None,
    frame_type: str = "braced",
    has_transverse_loading: bool = True,
    compression_sign: str = "negative",
) -> Dict[str, Any]:
    """Convenience function for Chapter H checks returning dict."""
    return aisc_chapter_h_check(
        loaded_beam,
        l_u,
        f_u,
        Ky=Ky,
        Kz=Kz,
        unbraced_positions_y=unbraced_positions_y,
        unbraced_positions_z=unbraced_positions_z,
        Cmx=Cmx,
        Cmy=Cmy,
        frame_type=frame_type,
        has_transverse_loading=has_transverse_loading,
        compression_sign=compression_sign,
    ).info()
