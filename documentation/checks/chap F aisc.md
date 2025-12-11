
# Scope (Chapter F — Beams and Channels)

This chapter covers allowable stresses for beams and channels. A brief preface aligning with the spec:

- Beams vs plate girders: treat a member as a “beam” for Chapter F when
   $$\frac{h}{t_w} \le \frac{970}{\sqrt{F_y}}$$
   where $h$ is the clear distance between flanges (in.), $t_w$ is the web thickness (in.), and $F_y$ is in ksi.
   If larger, bending stresses fall under Chapter G (plate girder provisions).
- Shear scope: use Chapter F for web shear checks unless designing for tension field action (see Chapter G).
- Applicability: applies to singly or doubly symmetric beams (including hybrids) and channels loaded in-plane
   through the shear center or braced against twist.
- Combined flexural + axial: see Section H1.

# F1 — Strong-Axis Bending (I-Shaped Members & Channels)


## F1.1 — Members with **Compact** Sections

**Use when**: I-shapes or channels, symmetric about and loaded in the plane of the minor axis; excludes hybrids and steels with $F_y>65$ ksi.

**Checklist**

* [ ] **Shape & loading**: I-shape or channel; symmetric; loading in plane of symmetry.
* [ ] **Compactness**: Section is **compact** per B5.1.
* [ ] **Material limits**: **Not hybrid**, and $F_y \le 65$ ksi.
* [ ] **Flange connection**: Flanges **continuously connected** to web(s).
* [ ] **Lateral bracing**: $L_b \le L_c$ (compression-flange unbraced length limit).

**Allowable stress & limit length**

* $F_b = 0.66\,F_y$.
* $L_c=\min\!\Big(\dfrac{76\,b_f}{\sqrt{F_y}},\ \dfrac{20{,}000}{(d/A_f)\,F_y}\Big)$. *(Units per spec; $F_y$ in ksi.)*

**Optional moment redistribution (where applicable)**

* Continuous over supports / rigidly framed: limited **$9/10$** reduction of **negative** gravity moments with compensating increase in **positive** moment; not for cantilevers. *(See spec paragraph for conditions.)*


## F1.2 — Members with **Noncompact** Sections

*(Still satisfy the F1.1 setup except compactness; three cases below.)*

### F1.2(a) — **Rolled** shapes, **noncompact flanges** (not built-up)

**Checklist**

* [ ] **Shape & loading**: I-shape or channel; symmetric; loading in plane of symmetry.
* [ ] **Compactness**: Flanges are **noncompact** per B5.1.
* [ ] **Material limits**: **Not hybrid**, and $F_y \le 65$ ksi.
* [ ] **Flange connection**: Flanges **continuously connected** to web(s).
* [ ] **Lateral bracing**: $L_b \le L_c$ (compression-flange unbraced length limit).
* [ ] **Rolled shape**: **Exclude** built-up members (see F1.2(b) for built-up).

**Allowable stress**

* $F_b = F_y\!\left[0.79 - 0.002\,\dfrac{b_f}{2t_f}\,\sqrt{F_y}\right]$.

---

### F1.2(b) — **Built-up** I-sections, **noncompact flanges**, web compact or noncompact

**Checklist**

* [ ] **Shape & loading**: Built-up I-shape; symmetric; loading in plane of symmetry.
* [ ] **Compactness**: Flanges are **noncompact** per B5.1; web is **compact or noncompact**.
* [ ] **Material limits**: **Not hybrid**, and $F_y \le 65$ ksi.
* [ ] **Flange connection**: Flanges **continuously connected** to web(s).
* [ ] **Lateral bracing**: $L_b \le L_c$ (compression-flange unbraced length limit).
* [ ] **Built-up member**: This path specifically applies to built-up (welded/bolted) sections.

**Allowable stress**

* $F_b = F_y\!\left[0.79 - 0.002\,\dfrac{b_f}{2t_f}\,\sqrt{F_y}\right]\dfrac{1}{k_c}$, with
  $k_c = 
  \begin{cases}
  1.0, & h/t_w \le 70\\[2pt]
  \dfrac{4.05}{(h/t_w)^{0.46}}, & h/t_w > 70
  \end{cases}$.  *(Per spec’s definition and bounds.)*

---

### F1.2(c) — **Other noncompact** sections not included above

**Checklist**

* [ ] Noncompact by B5, **not** covered in (a) or (b).
* [ ] **Loaded through the shear center**.
* [ ] **Laterally braced** in compression region at intervals $\le \dfrac{76\,b_f}{\sqrt{F_y}}$.

**Allowable stress**

* $F_b = 0.60\,F_y$.



## F1.3 — **Unbraced length $L_b > L_c$** (LTB governs)

> **Units per spec:** use **inches** for $l, r_T, d$; **in²** for $A_f$; **ksi** for $F_y$. $C_b \le 2.3$.

**Tension flange** (always, for $L_b > L_c$)

$$
F_{b,\;tension} = 0.60\,F_y \quad \text{(Eq. F1-5).}
$$

For **channels bent about their major axis**, take **compression** from **Eq. (F1-8) only**.

**Clarifications (per spec):**

* **Hybrid plate girders**: use **$F_y$ of the compression flange** in Eqs. **(F1-6)** and **(F1-7)**; **Eq. (F1-8) is not permitted**.
* **Tees**: Section **F1.3 does not apply** if the **stem is in compression** along the unbraced length.

---

### F1.3(a) — Inelastic LTB (Eq. F1-6)

**Use when the slenderness lies in this band**

$$
\sqrt{\frac{102\times 10^{3}\,C_b}{F_y}}
\;\le\;
\frac{l}{r_T}
\;\le\;
\sqrt{\frac{510\times 10^{3}\,C_b}{F_y}}.
$$

**Allowable compressive bending stress**

$$
F_{b,\;comp} \;=\;
\left[
\frac{2}{3}
\;-\;
\frac{F_y \,(l/r_T)^2}{1530\times 10^{3}\,C_b}
\right] F_y
\;\;\le\;\; 0.60\,F_y.
$$

**Checklist**

* [ ] $L_b > L_c$; I-shape/channel with axis of symmetry in, and loaded in, the **web plane**.
* [ ] $l, r_T, A_f, C_b$ computed (see “Extra” notes).
* [ ] Slenderness satisfies the **inelastic band** above.
* [ ] **Hybrid girders:** use **$F_y$ of the compression flange** in this formula; **do not** use Eq. (F1-8).

---

### F1.3(b) — Elastic LTB (Eq. F1-7)

**Use when slenderness exceeds the upper bound**

$$
\frac{l}{r_T} \;\ge\; \sqrt{\frac{510\times 10^{3}\,C_b}{F_y}}.
$$

**Allowable compressive bending stress**

$$
F_{b,\;comp} \;=\;
\frac{170\times 10^{3}\,C_b}{(l/r_T)^2}
\;\;\le\;\; 0.60\,F_y.
$$

**Checklist**

* [ ] $L_b > L_c$; symmetric I-shape/channel, loaded in web plane.
* [ ] $l, r_T, C_b$ computed; **$l/r_T$** above the elastic threshold.
* [ ] **Hybrid girders:** use **$F_y$ of the compression flange** here; **do not** use Eq. (F1-8).

---

### F1.3(c) — Local-flange/check (Eq. F1-8)

**Valid for any** $l/r_T$; **applies only if** the **compression flange is solid and approximately rectangular**, with **compression-flange area $\ge$** tension-flange area.
**Allowable compressive bending stress**

$$
F_{b,\;comp} \;=\;
\frac{12\times 10^{3}\,C_b}{\,l\,d/A_f\,}
\;\;\le\;\; 0.60\,F_y.
$$

**Checklist**

* [ ] Compression flange solid, \~rectangular; $A_{\text{comp}}\ge A_{\text{tension}}$.
* [ ] Compute with the same $l, d, A_f, C_b$ units noted above.
* [ ] **Hybrid girders:** **Eq. (F1-8) not permitted**.
* [ ] **Channels (major-axis bending):** use **Eq. (F1-8) only**.

---

### How to select the compression allowable (summary)

1. Compute $l, r_T, A_f, d, C_b$ (cap $C_b$ at 2.3; take $C_b=1.0$ if interior moment exceeds both end moments, or conservatively for cantilevers).
2. Check the **slenderness bounds** above:

   * If in the **inelastic band**, get $F_{b,\;comp}$ from **F1-6**.
   * If **beyond** the upper bound, get $F_{b,\;comp}$ from **F1-7**.
3. Also check **F1-8** (if permitted); **use the larger** of $\{\text{F1-6 or F1-7}\}$ and **F1-8**, but do not exceed $0.60F_y$ (and respect Chapter G where applicable).

**Definitions (per spec)**
$l$ = clear distance between brace points that **prevent twist and lateral displacement** of compression flange (cantilever braced only at the support may conservatively use the full cantilever length);
$r_T$ = radius of gyration of **compression flange + one-third of compression web** about an axis **in the web plane**;
$A_f$ = compression-flange area;
$C_b = 1.75 + 1.05(M_1/M_2) + 0.30(M_1/M_2)^2 \le 2.3$\*\*, with $M_1/M_2>0$ for reverse curvature, $<0$ for single curvature; use $C_b=1.0$ if the **maximum interior moment** in the unbraced length exceeds both end moments.



# F2 — Weak-Axis Bending (I-Shapes, Solid Bars & Rectangular Plates)

> Source: AISC (ASD) Chapter F, Section F2.

## F2.1 — Members with **Compact** Sections

**Use when**:

* Doubly-symmetric I-shapes with **compact flanges** (Sect. B5), bent about the weak axis.
* **Solid rounds, squares, or rectangles** bent about their weak axis.
* Flanges continuously connected to webs; excludes steels with $F_y>65$ ksi.

**Checklist**

* [ ] Shape: I-shape (doubly sym.), solid bar, or solid rectangle.
* [ ] **Compact** per B5.
* [ ] Flanges continuously connected to web.
* [ ] $F_y \le 65$ ksi.
* [ ] Bracing convenience note: Lateral bracing is not required for members loaded through the shear center about their weak axis,
      nor for members of equal strength about both axes.

**Allowable stress**

* $F_b = 0.75\,F_y$.



## F2.2 — Members with **Noncompact** Sections
### F2.2(a) — General noncompact case
**Checklist**

* [ ] **Shape**: Not box/tube/rectangular tube (see F3 for those); not I-beam with noncompact flanges (see F2.2(b)).
* [ ] **Compactness**: Section is **noncompact** per Sect. B5.
* [ ] **Flange connection**: Flanges continuously connected to web (for I-shapes).
* [ ] **Material limits**: $F_y \le 65$ ksi.
* [ ] **Loading**: Weak-axis bending.

**Allowable stress**

* **General case**: $F_b = 0.60\,F_y$.

### F2.2(b) — Doubly-symmetric I-shapes with noncompact flanges
**Checklist**

* [ ] **Shape**: Doubly-symmetric I-shape or H-shape.
* [ ] **Compactness**: Flanges are **noncompact** per Sect. B5.
* [ ] **Flange connection**: Flanges continuously connected to web.
* [ ] **Material limits**: $F_y \le 65$ ksi.
* [ ] **Loading**: Weak-axis bending.

**Allowable stress**

 $$
 F_b = F_y\!\left[1.075 - 0.005\left(\tfrac{b_f}{2t_f}\right)\sqrt{F_y}\right]
 $$



# F3 — Bending of Box Members, Rectangular Tubes & Circular Tubes

> Source: AISC (ASD) Chapter F, Section F3.

## F3.1 — Members with **Compact** Sections

**Use when**: Box/rectangular tube or circular tube; compact per Sect. B5.

**Checklist**

* [ ] Box or tube shape (rectangular or circular).
* [ ] Section **compact** (B5).
* [ ] Flanges continuously connected to webs.
* [ ] Meets **extra compactness rules**:

  * Depth $\le 6 \times$ width.
  * Flange thickness $\le 2 \times$ web thickness.
  * $L_b \le L_c$, with

    $$
    L_c = \frac{1950 + 1200(M_1/M_2)}{\sqrt{F_y}}
   \quad \text{not less than } 1200\left(\frac{b}{\sqrt{F_y}}\right).
    $$

**Allowable stress**

* $F_b = 0.66\,F_y$.

---

## F3.2 — Members with **Noncompact** Sections

**Checklist**

* [ ] **Shape**: Box, rectangular tube, square tube, or circular tube.
* [ ] **Compactness**: Section is **noncompact** by Sect. B5.
* [ ] **Flange connection**: Flanges continuously connected to webs.
* [ ] **Bracing check**:
  * If **depth < 6 × width**, lateral bracing not required.
  * If **depth ≥ 6 × width**, lateral support must be established by analysis.

**Allowable stress**

* $F_b = 0.60\,F_y$.

# F4 — Allowable Shear Stress

> **Units per spec:** use **inches** for $h,\ t_w,\ a$; **ksi** for $F_y$. This section sets the **allowable web shear stress** and when **stiffeners** are needed.

## Quick decision

1. **Stocky web**
   If $\dfrac{h}{t_w} \le \dfrac{380}{\sqrt{F_y}}$:

$$
F_v = 0.40\,F_y \quad \text{(on overall depth \(\times\) web thickness).}
$$

*(Eq. F4-1)*

2. **Slender web**
   If $\dfrac{h}{t_w} > \dfrac{380}{\sqrt{F_y}}$:

$$
F_v = C_v\,\big(0.40\,F_y\big) \quad \text{(on clear web depth \(h\) \(\times\) \(t_w\)).}
$$

*(Eq. F4-2; reduction via $C_v$)*

### $C_v$ and $k_v$ definitions (for F4)

The spec expresses $C_v$ using the following piecewise form in terms of slenderness and panel buckling:

$$
C_v =
\begin{cases}
\dfrac{45{,}000}{\dfrac{h}{t_w}\,\sqrt{k_v\,F_y}}, & \text{if } C_v \le 0.8 \\
\dfrac{190}{\left(\dfrac{h}{t_w}\right)\sqrt{F_y}}, & \text{if } C_v > 0.8
\end{cases}
$$

with the web-buckling coefficient $k_v$ defined by

$$
k_v =
\begin{cases}
4.00 + 5.34\,(a/h)^2, & a/h < 1.0 \\
5.34 + \dfrac{4.00}{(a/h)^2}, & a/h \ge 1.0
\end{cases}
$$

where $h$ = clear distance between flanges (in.), $t_w$ = web thickness (in.), and $a$ = clear distance between transverse stiffeners (in.).

## Variables (per spec)

* $t_w$: web thickness, **in.**
* $a$: clear distance between **transverse stiffeners**, **in.**
* $h$: clear distance between **flanges** at the section checked, **in.**

## Stiffeners & limits

* If $h/t_w>260$ **and** the computed web shear exceeds Eq. **F4-2**, **intermediate stiffeners are required** (see **F5**).
* **Maximum** $h/t_w$ limits and **alternative plate-girder shear (tension-field action)** are in **Chapter G**. For **coped ends**, check **J4**.

*(Source: AISC ASD, Chapter F.)*




# Extra

## Calculating $r_T$

**What $r_T$ is (spec):** radius of gyration of a section consisting of the **compression flange** plus **one-third of the compression web area**, taken **about an axis in the plane of the web** (i.e., the vertical axis through the web mid-plane).

### Short, accurate method for W-shapes (strong-axis bending)

1. **Clear web depth**

   $$
   h \;\approx\; d - 2t_f 
   \quad\text{(better: use the clear distance between flange fillets).}
   $$

2. **Build the “compression tee”**
   Flange area:

   $$
   A_f = b\,t_f
   $$

   One-third of the compression web area:

   $$
   A_{w,c}=\frac{t_w\,h}{3}\cdot\frac{1}{3}=\frac{t_w\,h}{6}
   $$

   Tee area:

   $$
   A_T = A_f + A_{w,c} = b t_f + \frac{t_w h}{6}
   $$

3. **Second moments about the vertical axis (in the web plane)**
   (Both centroids lie on this axis, so no parallel-axis terms are needed.)

   $$
   I_{y,f}=\frac{t_f\,b^{3}}{12}, 
   \qquad 
   I_{y,w}=\frac{(h/6)\,t_w^{3}}{12}=\frac{h\,t_w^{3}}{72}
   $$

4. **Compute $r_T$**

   $$
   r_T
   \;=\;
   \sqrt{\frac{I_{y,f}+I_{y,w}}{A_T}}
   \;=\;
   \sqrt{\frac{\dfrac{t_f b^{3}}{12}+\dfrac{h\,t_w^{3}}{72}}
                {\,b t_f+\dfrac{t_w h}{6}\,}}
   $$

5. **Useful approximation (usually tight for rolled W-shapes)**

   $$
   r_T \;\approx\;
   \sqrt{\frac{t_f b^{3}/12}{b t_f}}
   \;=\;\frac{b}{\sqrt{12}}
   \;\approx\;0.289\,b
   \quad\text{(the web’s }t_w^{3}\text{ term is tiny).}
   $$

> **Spec alignment:** This matches the AISC definition of $r_T$ (compression flange $+$ one-third compression web, about an axis **in the web plane**). For built-up or unusual geometry, compute the exact compression-tee section properties from dimensions/fillets rather than using $h\approx d-2t_f$.

---

## Calculating $C_b$

**What $C_b$ is (spec):** the **moment gradient factor** that amplifies (or reduces) the lateral-torsional buckling capacity for **unbraced lengths**; bounded above by **2.3**. Use **$C_b=1.0$** when the **maximum moment within the unbraced length exceeds the end moments**, when frames are **braced against joint translation** for $F_{bx}$ checks, and it is always conservative to take $C_b=1.0$. For **cantilevers**, $C_b=1.0$ is a conservative default.

### Formula and sign convention

$$
C_b \;=\; 1.75 \;+\; 1.05\left(\frac{M_1}{M_2}\right)
           \;+\; 0.30\left(\frac{M_1}{M_2}\right)^2,
\qquad C_b \le 2.3,
$$

where $M_1$ is the **smaller** and $M_2$ the **larger** **end moment** (by magnitude) over the unbraced segment, about the **strong axis**.

* **Sign of $M_1/M_2$:** **positive** when end moments have the **same sign** (*reverse-curvature* bending in the spec’s terminology) and **negative** when they have **opposite signs** (*single-curvature*).
* If any **interior point moment** exceeds both end moments, **take $C_b=1.0$**.
* **Clip** any computed $C_b$ to **2.3**.

### Quick reference (common cases)

* **Uniform moment** (equal end moments, same sign): $M_1/M_2=+1\Rightarrow C_b=3.10 \rightarrow \text{use }2.3$.
* **Simply supported w/ midspan maximum** (e.g., UDL or point load at midspan): interior moment $>$ ends $\Rightarrow C_b=1.0$.
* **Equal and opposite end moments** (pure single-curvature): $M_1/M_2=-1 \Rightarrow C_b=1.0$.
* **Cantilever**: conservative **$C_b=1.0$**.

> **Spec alignment:** The formula, sign convention, special cases (interior maximum, frames braced against joint translation, cantilevers), and cap at 2.3 are straight from AISC Chapter F1.3 notes. Use the same $C_b$ value consistently with the chosen LTB equation (F1-6 or F1-7).

---

**Tip:** In your F1.3 checklist, group the LTB inputs as: $l$, $r_T$ (from above), $A_f$, and $C_b$ (from above), then evaluate **(F1-6) or (F1-7)** and **(F1-8)** as applicable; take the **larger** allowable compressive stress, observing any caps in Chapter G.
