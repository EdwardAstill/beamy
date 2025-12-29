

# Chapter H — Combined Stresses (ASD)

> Source: **AISC ASD Specification, Chapter H — Combined Stresses** 

This chapter governs **members subjected to combined axial force and bending**.

**Scope**

* Applies to **doubly and singly symmetric members**
* Uses:

  * **Chapter E** for allowable axial stress (F_a)
  * **Chapter F** for allowable bending stresses (F_{bx}, F_{by})
* Covers:

  * Axial **compression + bending**
  * Axial **tension + bending**

---

## H1 — Axial Compression and Bending

Members subjected to **axial compression** and **bending about one or both axes** must satisfy interaction requirements at all critical sections.

### Stress Notation

* (f_a) = computed axial compressive stress
* (f_{bx}, f_{by}) = computed compressive bending stresses
* (F_a) = allowable axial compressive stress (Chapter E)
* (F_{bx}, F_{by}) = allowable bending stresses (Chapter F)
* (F'*{ex}, F'*{ey}) = Euler buckling stress (with ASD safety factor)
* (C_{mx}, C_{my}) = moment gradient coefficients

Subscripts:

* (x, y) → bending axis
* (b) → bending
* (a) → axial
* (e) → Euler buckling

---

### H1.1 — General Interaction Equations

Members must satisfy **both** of the following:

#### (H1-1) Buckling-Modified Interaction

[
\frac{f_a}{F_a}
+
\frac{C_{mx} f_{bx}}{\left(1 - \dfrac{f_a}{F'*{ex}}\right) F*{bx}}
+
\frac{C_{my} f_{by}}{\left(1 - \dfrac{f_a}{F'*{ey}}\right) F*{by}}
;\le; 1.0
]

#### (H1-2) Linear Interaction Check

[
\frac{f_a}{0.60F_y}
+
\frac{f_{bx}}{F_{bx}}
+
\frac{f_{by}}{F_{by}}
;\le; 1.0
]

> **Both equations must be satisfied** unless the simplified provision below applies.

---

### H1.2 — Simplified Interaction (Low Axial Load)

When:

[
\frac{f_a}{F_a} ;\le; 0.15
]

Equation **(H1-3)** may be used **instead of** (H1-1) and (H1-2):

[
\frac{f_a}{F_a}
+
\frac{f_{bx}}{F_{bx}}
+
\frac{f_{by}}{F_{by}}
;\le; 1.0
]

**Notes**

* No moment-gradient modifiers
* No Euler amplification terms
* Valid only for **small axial compression**

---

### H1.3 — Euler Stress for Interaction

[
F'_e
====

\frac{12\pi^2 E}{23,(K l_b / r_b)^2}
]

where:

* (l_b) = unbraced length **in the plane of bending**
* (r_b) = radius of gyration for that plane
* (K) = effective length factor for that plane

> As with (F_a) and (F_b), (F'_e) may be increased by **1/3** per Section A5.2.

---

### H1.4 — Moment Gradient Coefficient (C_m)

The coefficient (C_m) accounts for moment shape along the unbraced length.

#### (a) Frames **with sidesway**

[
C_m = 0.85
]

---

#### (b) Braced frames, **no transverse loading**, rotationally restrained

[
C_m = 0.6 - 0.4\left(\frac{M_1}{M_2}\right)
]

where:

* (M_1) = smaller end moment
* (M_2) = larger end moment
* (M_1/M_2 > 0): reverse curvature
* (M_1/M_2 < 0): single curvature

---

#### (c) Braced frames **with transverse loading**

* (C_m) may be obtained from analysis
* If not analyzed, use:

| End restraint condition            | (C_m) |
| ---------------------------------- | ----- |
| Ends restrained against rotation   | 0.85  |
| Ends unrestrained against rotation | 1.0   |

---

## H2 — Axial Tension and Bending

Members subjected to **axial tension** and **bending** must satisfy:

[
\frac{f_a}{F_t}
+
\frac{f_{bx}}{F_{bx}}
+
\frac{f_{by}}{F_{by}}
;\le; 1.0
\quad \text{(H2-1)}
]

where:

* (f_a) = computed axial **tensile** stress
* (f_b) = computed **tensile** bending stress
* (F_t) = allowable tensile stress (Section D1)

### Important Limitation

Any **compressive bending stress** caused by a *different* load source (relative to the axial tension) must **also** satisfy the applicable limits of **Chapter F**.

---

## Practical Design Workflow (Chapter H)

1. Determine:

   * (f_a) from axial force
   * (f_{bx}, f_{by}) from bending moments
2. Obtain:

   * (F_a) from Chapter E
   * (F_{bx}, F_{by}) from Chapter F
3. Compute:

   * (F'*{ex}, F'*{ey})
   * (C_{mx}, C_{my})
4. Check:

   * (H1-1) **and** (H1-2), or
   * (H1-3) if (f_a/F_a \le 0.15)
5. For tension members, use (H2-1)

---

## Key Observations

* Chapter H is **stress-based**, not force-based
* Buckling interaction is handled explicitly via (F'_e)
* Moment gradient effects matter in compression
* Simplified interaction is intentionally conservative and limited

---

If you want, I can:

* Convert this to a `.md` file matching your repo style
* Add **worked numerical examples**
* Map these checks directly to your **beamy** implementation
* Provide a **flowchart** version for quick hand checks
