C2. FRAME STABILITY 
1. Braced Frames 
In trusses and in those frames where lateral stability is provided by adequate 
attachment to diagonal bracing, to shear walls, to an adjacent structure having 
adequate lateral stability or to floor slabs or roof decks secured horizontally by 
walls or bracing systems parallel to the plane of the frame, the effective length 
factor K for the compression members shall be taken as unity, unless analysis 
shows that a smaller value is permitted. 

2. Unbraced Frames 
In frames where lateral stability is dependent upon the bending stiffness of rig-
idly connected beams and columns, the effective length Kl of compression 
members shall be determined by analysis and shall not be less than the actual 
unbraced length.

# Chapter E — Columns and Other Compression Members (ASD)

> Source: AISC 360 (ASD), Chapter E — Columns and Other Compression Members.

This chapter applies to **prismatic compression members** with **compact or noncompact** cross-sections,
subject to **axial compression through the centroidal axis**.

**Not covered here:**
- Members with **slender elements** → see **Appendix B5.2**
- Members with **combined axial compression and flexure** → see **Chapter H**
- **Tapered members** → see **Appendix F7**

---

# E1 — Effective Length and Slenderness Ratio

## Applicability
- Axially loaded compression members
- Compression through centroidal axis
- Governing buckling mode is **flexural**, unless noted otherwise

## Definitions

**Effective length**
\[
L_e = K L
\]

where:
- \(L\) = unbraced length
- \(K\) = effective-length factor (per Section C2)

**Slenderness ratio**
\[
\frac{K L}{r}
\]

where:
- \(r\) = radius of gyration about the governing buckling axis

### Governing slenderness
- Use the **largest** \(KL/r\) of any unbraced segment
- Use the **least** radius of gyration

### Limiting slenderness
- Maximum permitted \(KL/r\) limits are given in **Section B7**

---

# E2 — Allowable Axial Compressive Stress

## Elastic–Inelastic Transition Parameter

\[
C_c = \sqrt{\frac{2\pi^2 E}{F_y}}
\]

where:
- \(E = 29{,}000\) ksi
- \(F_y\) = yield stress (ksi)

---

## E2.1 — Inelastic Buckling Range  
### \( \dfrac{K L}{r} < C_c \)

**Allowable compressive stress**

\[
F_a
=
\frac{
\left[
1 - \dfrac{(KL/r)^2}{2C_c^2}
\right] F_y
}{
\frac{5}{3}
+
\frac{3(KL/r)}{8C_c}
-
\frac{(KL/r)^3}{8C_c^3}
}
\quad \text{(Eq. E2-1)}
\]

### Checklist
- [ ] Axial compression only
- [ ] Cross-section meets **Table B5.1**
- [ ] \(KL/r < C_c\)
- [ ] Use **gross area**

---

## E2.2 — Elastic Buckling Range  
### \( \dfrac{K L}{r} \ge C_c \)

**Allowable compressive stress**

\[
F_a
=
\frac{12\pi^2 E}{23 (KL/r)^2}
\quad \text{(Eq. E2-2)}
\]

### Checklist
- [ ] Axial compression only
- [ ] Governing slenderness ratio exceeds \(C_c\)
- [ ] Euler-type buckling controls

---

## Summary — Selecting \(F_a\)

1. Compute \(KL/r\) for each principal axis
2. Determine governing (largest) value
3. Compute \(C_c\)
4. Use:
   - **Eq. E2-1** if \(KL/r < C_c\)
   - **Eq. E2-2** if \(KL/r \ge C_c\)

---

# E3 — Flexural–Torsional Buckling

## When required

Flexural–torsional or torsional buckling must be considered for:

- Singly symmetric sections:
  - Angles
  - Tees
- Unsymmetric sections
- Certain built-up sections with **thin elements**

### Notes
- Doubly symmetric rolled shapes usually buckle flexurally
- Open sections are more susceptible to torsional instability
- See **Appendix E3** or perform special analysis when applicable

---

# E4 — Built-Up Compression Members

## General Requirements

- All components must satisfy **Section B7**
- Fastener spacing must permit **force transfer**
- All end components in contact must be:
  - Fully bolted (tightened per Table J3.7), or
  - Continuously welded

---

## E4.1 — Fastener Spacing

### Maximum longitudinal spacing
- Rolled shapes in contact: **24 in**
- Plate components:
  - Unstaggered:  
    \[
    \le \min\left(127/\sqrt{F_y}\;t,\;12\text{ in}\right)
    \]
  - Staggered:
    \[
    \le \min\left(190/\sqrt{F_y}\;t,\;18\text{ in}\right)
    \]

---

## E4.2 — Component Slenderness

- Slenderness of any component between connectors:
\[
\left(\frac{KL}{r}\right)_{\text{component}}
\le
\frac{3}{4}
\left(\frac{KL}{r}\right)_{\text{member}}
\]

- Use **least** component radius of gyration
- Minimum **two intermediate connectors** required

---

## E4.3 — Lacing and Tie Plates

### Lacing requirements
- Resist transverse shear = **2%** of axial force
- Slenderness limits:
  - Single lacing: \(L/r \le 140\)
  - Double lacing: \(L/r \le 200\)

### Geometry
- Inclination:
  - ≥ 60° (single)
  - ≥ 45° (double)

### Tie plates
- Required at:
  - Both ends
  - Any interruption in lacing
- Thickness ≥ \(1/50\) of fastener-line spacing

---

# E5 — Pin-Connected Compression Members

- Must comply with **Section D3**
- Member assumed pinned for buckling analysis
- Connection detailing must ensure free rotation

---

# E6 — Column Web Shear

- Column webs must be checked for **concentrated force introduction**
- Use **Section K1** for:
  - Beam reactions
  - Bracket loads
  - Girder framing forces

---

# Design Workflow (Practical)

1. Determine **unbraced lengths**
2. Select **effective length factor \(K\)** (or use Direct Analysis)
3. Compute governing \(KL/r\)
4. Compute \(C_c\)
5. Determine allowable \(F_a\)
6. Check:
\[
\frac{P}{A_g} \le F_a
\]

---

# Notes for FEM / Frame Analysis

- \(K\) reflects **end rotational restraint**, not just fixity flags
- Spring restraints can be used to estimate effective fixity
- **Direct Analysis Method** eliminates explicit \(K\)-factors
- Extract **axial force demand** from global model, then check member per Chapter E

---

# beamy implementation

The axial compression check is implemented in the codebase as:

- [src/beamy/checks/aisc_9.py](src/beamy/checks/aisc_9.py)

Usage on a 1D extracted member (a `LoadedMember`):

- `loaded_beam.check_aisc_chapter_e(length_unit="m", force_unit="N", K=1.0)`
- Optional bracing segmentation: `unbraced_positions=[0.0, 0.9, 1.8]`

This check is **axial-only** (no P–M interaction). For combined axial + flexure, use Chapter H.

---
