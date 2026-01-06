
## 1️⃣ Releases

> **“Releases is about how a member is connected to another – it is a member attribute.”**

✔ **Correct**, with one precision tweak:

> **Releases describe how a *member end transmits stiffness* at a connection.**

Key points:

* A release does **not** decide *whether* a connection exists.
* It decides *what force components that member can carry* across that connection.
* It applies **per member end**, not per joint.

So this is right:

* member-level
* end-specific
* stiffness-related

### Implementation: Static Condensation

For a professional FEA implementation, use **Static Condensation** rather than soft springs.

1.  **Element Level**: Perform condensation inside `beam3d.py`.
2.  **Partition** the element stiffness matrix $K_e$ into retained ($r$) and released ($c$) degrees of freedom.
3.  **Condense** the matrix:
    $$K_{condensed} = K_{rr} - K_{rc} K_{cc}^{-1} K_{cr}$$
4.  **Store** the condensation matrix ($K_{cc}^{-1} K_{cr}$) on the Element for instant recovery.
5.  **Assemble** $K_{condensed}$ into the global system.
6.  **Recover** displacements for the released DOFs using the stored matrix after solving global $u$.

This ensures the released behavior is mathematically exact and avoids numerical instability.

---

## 2️⃣ Connections

> **“Connections is about which members are connected to each other – it is a frame attribute.”**

✔ **Exactly correct.**

Why:

* Connections reference **multiple entities**
* They describe **topology / compatibility**
* They answer: *who shares motion with whom?*

Examples:

* two member ends share a node
* a member end connects to the middle of another member
* rigid link / MPC between two nodes

All of those are:

* frame-level
* independent of stiffness behavior

---

## 3️⃣ Supports

> **“Supports is about how a node is supported – it is a load case attribute.”**

✔ **Correct** — and this is the one people most often get wrong.

Why supports must be load-case-level:

* supports can change between load cases
* construction staging exists
* settlement / imposed displacement exists
* nonlinear analysis often needs different supports per step

Supports answer:

> “Which DOFs are restrained to ground *for this analysis case*?”

They are **boundary conditions**, not geometry or element properties.

---

## The clean final mental model (lock this in)

You can summarize the entire architecture like this:

| Concept        | Question it answers                                  | Lives on   |
| -------------- | ---------------------------------------------------- | ---------- |
| **Connection** | Who shares motion with whom?                         | `Frame`    |
| **Release**    | What stiffness does this member transmit at its end? | `Member`   |
| **Support**    | Which DOFs are restrained to ground for this case?   | `LoadCase` |

None of these imply the others.

**Implementation note (aligns with `idea.md`)**

- **DOF order (3D)**: `[UX, UY, UZ, RX, RY, RZ]`
- **Release masks** (on `Member`): True = released (that member end does *not* transmit that component)
- **Restraint masks** (on `Support`): True = restrained (fixed to ground for that load case)

---

## One last clarification (important for second-order)

A release is **not** “how a member is connected to another” in a geometric sense.
It is **how a member participates in the joint equilibrium**.

Two members can:

* be fully connected
* share the same node
* and yet one of them transmits *no moment* at that joint

That’s the whole point of releases.

---

## Verdict on your model

Your architecture is now:

* **mechanically sound**
* **scales to second-order**
* **matches professional FEA practice**
* **cleanly separates concerns**

You’re no longer fighting the physics with the API — which is the hard part.

If you want, next we can:

* map this directly to the global stiffness assembly step,
* define the minimal 2D DOF version (much simpler),
* or design the validation rules that make the API feel “smart” without guessing physics.

Just say the word.
