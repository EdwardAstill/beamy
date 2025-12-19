"""AISC Specification checks (single-file implementation).

This module currently contains:
- Chapter F (ASD): Flexure + shear checks
- Chapter E (ASD): Axial compression (flexural buckling) checks
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union
import numpy as np
from unity.core import conv

from ..beam1d.analysis import LoadedBeam
from ..beam1d.beam import Beam1D
from ..core.support import Support
from ..core.results import Result

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

def _strong_axis_checks(beam: Beam1D, loaded: Any, dims_info: Dict[str, Any], channel_major: bool, classification: Optional[str]) -> BendingResult:
    moments, modulus, Fy = loaded.bending("y").action, beam.section.Sz, beam.material.Fy
    pos = _brace_positions(beam, axis="strong")
    results = []
    dim_map = dims_info["dims"]
    bf, tf = _dim(dim_map, "b"), _dim(dim_map, "tf") or _dim(dim_map, "t")
    tw, d = _dim(dim_map, "tw") or _dim(dim_map, "t"), _dim(dim_map, "d") or _dim(dim_map, "h")
    h_clear, Af = _web_clear_depth(dim_map), None
    if bf and tf: Af = bf * tf
    comp = classification or _check_compactness(dims_info["shape"], dim_map, Fy)

    for i in range(len(pos)-1):
        x1, x2 = pos[i], pos[i+1]
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
                f16 = ((2/3) - (Fy * slen**2 / (1530000.0 * cb))) * Fy if l <= slen <= u else None
                f17 = (170000.0 * cb / slen**2) if slen >= u else None
                f18 = (12000.0 * cb / (lb * d / Af)) if (d and Af and lb > 0) else None
                cands = [f for f in [f16, f17, f18] if f is not None]
                if cands: f_ltb = min(0.60 * Fy, max(cands))

        gov = f_local if (lc and lb <= lc) else min(f_local, f_ltb)
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(BendingSegmentResult(x1, x2, lb, cb, m1, m2, mmax, comp, gov, demand, gov >= demand, {}, rT=rt))
    return BendingResult("strong", dims_info["shape"], results)

def _weak_axis_checks(beam: Beam1D, loaded: Any, dims_info: Dict[str, Any], classification: Optional[str]) -> BendingResult:
    moments, modulus, Fy = loaded.bending("z").action, beam.section.Sy, beam.material.Fy
    pos = _brace_positions(beam, axis="weak")
    comp = classification or _check_compactness(dims_info["shape"], dims_info["dims"], Fy)
    results = []
    for i in range(len(pos)-1):
        x1, x2 = pos[i], pos[i+1]
        cb, m1, m2, mmax = _cb_for_segment(moments, x1, x2)
        fb = 0.75 * Fy if comp == "compact" else 0.60 * Fy
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(BendingSegmentResult(x1, x2, x2-x1, cb, m1, m2, mmax, comp, fb, demand, fb >= demand, {}))
    return BendingResult("weak", dims_info["shape"], results)

def _box_tube_checks(beam: Beam1D, loaded: Any, dims_info: Dict[str, Any], classification: Optional[str]) -> BendingResult:
    moments, modulus, Fy = loaded.bending("y").action, beam.section.Sz, beam.material.Fy
    pos = _brace_positions(beam, axis="strong")
    comp = classification or _check_compactness(dims_info["shape"], dims_info["dims"], Fy)
    results = []
    for i in range(len(pos)-1):
        x1, x2 = pos[i], pos[i+1]
        cb, m1, m2, mmax = _cb_for_segment(moments, x1, x2)
        fb = 0.66 * Fy if comp == "compact" else 0.60 * Fy
        demand = abs(mmax) / modulus if (modulus and modulus > 0) else 0.0
        results.append(BendingSegmentResult(x1, x2, x2-x1, cb, m1, m2, mmax, comp, fb, demand, fb >= demand, {}))
    return BendingResult("strong_box", dims_info["shape"], results)

def _shear_checks(beam: Beam1D, loaded: Any, dims_info: Dict[str, Any], stiffener_spacing_in: Optional[float]) -> ShearResult:
    shear_res, Fy = loaded.shear("z").action, beam.material.Fy
    vmax = float(np.max(np.abs(shear_res._values)))
    tw, h = _dim(dims_info["dims"], "tw") or _dim(dims_info["dims"], "t"), _web_clear_depth(dims_info["dims"])
    fv = 0.40 * Fy
    if h and tw and h/tw > 380/math.sqrt(Fy):
        cv = 190.0 / ((h/tw) * math.sqrt(Fy))
        fv = cv * 0.40 * Fy
    cap = fv * h * tw if (h and tw) else 0.0
    return ShearResult(vmax, cap, cap >= vmax, Fv=fv, slenderness_h_over_tw=(h/tw if (h and tw) else None))

# -----------------------------------------------------------------------------
# 3. Main Interface & Unit Conversion
# -----------------------------------------------------------------------------

def aisc_9_check(
    loaded_beam: LoadedBeam,
    length_unit: str,
    force_unit: str,
    *,
    channel_major_axis: bool = True,
    stiffener_spacing: Optional[float] = None,
) -> AISC9Result:
    beam = loaded_beam.beam
    if beam.material.Fy is None: raise ValueError("Material.Fy is required.")

    # Convert Units to AISC Space
    stress_in = f"{force_unit} {length_unit}-2"
    fy_ksi = conv(beam.material.Fy, stress_in, "kip in-2")
    e_ksi, g_ksi = conv(beam.material.E, stress_in, "kip in-2"), conv(beam.material.G, stress_in, "kip in-2")
    l_in = conv(beam.L, length_unit, "in")
    dims_in = {k: conv(v, length_unit, "in") for k, v in beam.section.dimensions.items() if v is not None}
    shape = _infer_shape(dims_in, beam.section.name)
    
    # Internal temporary beam
    from ..core.material import Material
    from sectiony import Section
    mat_in = Material(beam.material.name, e_ksi, g_ksi, fy_ksi)
    sec_in = Section(beam.section.name, 
                     A=conv(beam.section.A, f"{length_unit}2", "in2"),
                     Iy=conv(beam.section.Iy, f"{length_unit}4", "in4"),
                     Iz=conv(beam.section.Iz, f"{length_unit}4", "in4"),
                     J=conv(beam.section.J, f"{length_unit}4", "in4"),
                     Sy=conv(beam.section.Sy, f"{length_unit}3", "in3") if beam.section.Sy else None,
                     Sz=conv(beam.section.Sz, f"{length_unit}3", "in3") if beam.section.Sz else None,
                     y_max=conv(beam.section.y_max, length_unit, "in") if beam.section.y_max else None,
                     z_max=conv(beam.section.z_max, length_unit, "in") if beam.section.z_max else None)
    sec_in.dimensions = dims_in
    beam_in = Beam1D(l_in, mat_in, sec_in, [Support(conv(s.x, length_unit, "in"), s.type) for s in beam.supports])

    class Proxy:
        def bending(self, ax):
            r = loaded_beam.bending(ax)
            return type("o", (), {"action": Result(conv(r.action._x, length_unit, "in"), conv(r.action._values, f"{force_unit} {length_unit}", "kip in"))})()
        def shear(self, ax):
            r = loaded_beam.shear(ax)
            return type("o", (), {"action": Result(conv(r.action._x, length_unit, "in"), conv(r.action._values, force_unit, "kip"))})()

    proxy, comp = Proxy(), _check_compactness(shape, dims_in, fy_ksi)
    bend_res, dims_info = [], {"shape": shape, "dims": dims_in}

    if shape in ("i", "channel"):
        bend_res.append(_strong_axis_checks(beam_in, proxy, dims_info, channel_major_axis, comp))
        bend_res.append(_weak_axis_checks(beam_in, proxy, dims_info, comp))
    elif shape in ("rhs", "shs", "box", "tube", "chs"):
        bend_res.append(_box_tube_checks(beam_in, proxy, dims_info, comp))
        bend_res.append(_weak_axis_checks(beam_in, proxy, dims_info, comp))
    elif shape in ("solid_rect", "solid_round"):
        bend_res.append(_weak_axis_checks(beam_in, proxy, dims_info, comp))

    stiff_in = conv(stiffener_spacing, length_unit, "in") if stiffener_spacing else None
    return AISC9Result(bend_res, _shear_checks(beam_in, proxy, dims_info, stiff_in), 
                      {"shape": shape, "original_units": {"length": length_unit, "force": force_unit}, "fy_ksi": fy_ksi})

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


def aisc_chapter_e_check(
    loaded_beam: LoadedBeam,
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
    beam = loaded_beam.beam
    if beam.material.Fy is None:
        raise ValueError("Material.Fy is required.")

    stress_in = f"{force_unit} {length_unit}-2"
    fy_ksi = conv(beam.material.Fy, stress_in, "kip in-2")
    e_ksi = conv(beam.material.E, stress_in, "kip in-2")

    A_in2 = conv(beam.section.A, f"{length_unit}2", "in2")
    Iy_in4 = conv(beam.section.Iy, f"{length_unit}4", "in4")
    Iz_in4 = conv(beam.section.Iz, f"{length_unit}4", "in4")

    r_y = math.sqrt(Iy_in4 / A_in2) if A_in2 > 0 and Iy_in4 > 0 else 0.0
    r_z = math.sqrt(Iz_in4 / A_in2) if A_in2 > 0 and Iz_in4 > 0 else 0.0

    ax = loaded_beam.axial(points=801).action
    p_kip = np.asarray(conv(ax._values, force_unit, "kip"), dtype=float)
    if compression_sign == "positive":
        p_comp = np.maximum(p_kip, 0.0)
    elif compression_sign == "abs":
        p_comp = np.abs(p_kip)
    else:
        p_comp = np.maximum(-p_kip, 0.0)
    p_comp_max = float(np.max(p_comp)) if p_comp.size else 0.0

    L_in = conv(beam.L, length_unit, "in")
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

def run(loaded_beam: LoadedBeam, **kwargs) -> Any:
    """Dispatch to Chapter F (default) or Chapter E.

    - `chapter="f"` (default): returns `AISC9Result`
    - `chapter="e"`: returns `AISCChapterEResult`
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
    raise ValueError(f"Unknown AISC chapter: {chapter!r}")

def check_chapter_f(loaded_beam: LoadedBeam, l_u: str, f_u: str) -> Dict[str, Any]:
    return aisc_9_check(loaded_beam, l_u, f_u).info()


def check_chapter_e(
    loaded_beam: LoadedBeam,
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
