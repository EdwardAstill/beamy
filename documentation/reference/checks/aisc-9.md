## Overview
Summary of the AISC ASD 9th Edition Chapter F checks implemented in `src/beamy/checks/aisc_9.py`. All units are converted internally to ksi and inches.

## Inputs
- `loaded_beam` (analysis result): solved `LoadedBeam` with bending and shear actions.
- Units: `length_unit`, `force_unit` describing the model units.
- Optional switches:
  - `channel_major_axis` (default True): channels with Lb > Lc use Eq. F1-8 only.
  - `stiffener_spacing`: clear spacing for shear kv (same units as model). If omitted, kv defaults to 5.34 (unstiffened panel).

## Simplifications/Assumptions
- All members are treated as **rolled**; built-up provisions and kc factors are not considered.
- Chapter F scope only; combined axial/flexure and Chapter G plate-girder provisions are not implemented.
- Basic compact / noncompact / slender classification is performed using B5.1-style limits for flanges, webs, and box/tube/CHS walls.
- Weak-axis compactness for solid bars/rectangles is treated as compact by default.
- Slender elements are **not** treated with full AISC slender-element formulas. Instead, their bending allowable is limited to 0.60 Fy and a `scope_warning` is issued.
- Cb is capped at 2.3; Cb = 1.0 if interior moment exceeds both ends.
- Channels bent about major axis default to Eq. F1-8 only (may be relaxed via `channel_major_axis=False`).
- For box/tube compact checks (F3.1) on RHS/SHS/box: require h ≤ 6b, tf ≤ 2tw, and Lb ≤ Lc before using 0.66 Fy; otherwise 0.60 Fy is used.
- Circular hollow sections (CHS): compact when D/t ≤ 3300/Fy (ksi) and checked at 0.66 Fy; for D/t above this limit, treated as slender and checked at 0.60 Fy with a scope warning.
- Fy > 65 ksi and web slenderness h/tw > 970/√Fy are flagged in `scope_warnings` but not blocked.

## Checks Implemented
- **Strong-axis bending (F1)**: I-shapes and channels, with compact/noncompact/local flange/LTB equations; channels follow F1-8 only when `channel_major_axis` is True.
- **Weak-axis bending (F2)**: I-shapes, channels, solid rectangles/rounds; compact (0.75Fy) and noncompact paths.
- **Box/tube bending (F3)**: RHS/SHS/box/CHS.
  - RHS/SHS/box: compact 0.66 Fy when h ≤ 6b, tf ≤ 2tw, and Lb ≤ Lc; otherwise 0.60 Fy.
  - CHS: compact for D/t ≤ 3300/Fy and checked at 0.66 Fy; for D/t > 3300/Fy treated as slender and checked at 0.60 Fy with a scope warning.
- **Web shear (F4)**: Computes Fv, Cv, kv when spacing provided (or kv = 5.34 if spacing is omitted/unstiffened); reports capacity and stiffener need indicator (h/tw > 260 and capacity < Vmax).

## Output
```python
{
  "inputs": {
    "shape": <str>,
    "original_units": {
      "length": <str>,
      "force": <str>
    },
    "aisc_units": {
      "length": "in",
      "force": "kip",
      "stress": "ksi"
    },
    "material_Fy_ksi": <float>,
    "E_ksi": <float>,
    "G_ksi": <float>,
    "section_props_in": {
      "A_in2": <float>,
      "Iy_in4": <float>,
      "Iz_in4": <float>,
      "Sy_in3": <float | None>,
      "Sz_in3": <float | None>,
      "J_in4": <float>
    },
    "section_dims_in": {
      <dim_key>: <float>,
      ...
    },
    "scope_warnings": [<str>, ...]        # optional
  },

  "bending": [
    {
      "axis": <"strong" | "weak" | "strong_box">,
      "shape": <str>,
      "segments": [
        {
          "segment": {
            "x_start": <float>,
            "x_end": <float>,
            "Lb": <float>,
            "Cb": <float>
          },
          "moments": {
            "M1": <float>,
            "M2": <float>,
            "Mmax": <float>
          },
          "rT": <float | None>,
          "classification": <"compact" | "noncompact" | "slender">,
          "allowables": {
            "F1.1_compact": { "Fb": <float | None>, "Lc": <float | None> },
            "F1.2a_noncompact_rolled": { "Fb": <float | None> },
            "F1.2c_other": { "Fb": <float | None> },
            "F1.3a_inelastic": { "Fb": <float | None> },
            "F1.3b_elastic": { "Fb": <float | None> },
            "F1.3c_local_flange": { "Fb": <float | None> }
          },
          "demand": {
            "Fb_required": <float>
          },
          "pass": <bool>,
          "governing_Fb": <float>
        }
      ]
    }
  ],

  "shear": {
    "V_max": <float>,
    "slenderness_h_over_tw": <float | None>,
    "Cv": <float | None>,
    "kv": <float | None>,
    "Fv": <float | None>,
    "capacity": <float | None>,
    "a_over_h": <float | None>,
    "stiffeners_required": <bool | None>,
    "pass": <bool>
  }
}
```