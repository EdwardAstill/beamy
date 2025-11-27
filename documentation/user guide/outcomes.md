analyse should produce the folowing results a numpy array or list of tuples:

## Theory

This is how beams behave:

| Behaviour      | Action     | Stress     | Displacement |
| -------------- | ---------- | ---------- | ------------ |
| **Transverse** | V(x), M(x) | τ(x), σ(x) | **w(x)**     |
| **Axial**      | N(x)       | σ(x)       | u(x)         |
| **Torsion**    | T(x)       | τ(x)       | θ(x)         |

## Main analysis implementation
For my analysis i will break it down liike so:
| Behaviour      | Action     | Stress     | Displacement |
| -------------- | ---------- | ---------- | ------------ |
| **Shear**      | V(x)       | τ(x)       |              |
| **Bending**    | M(x)       | σ(x)       |              |
| **Axial**      | N(x)       | σ(x)       | u(x)         |
| **Torsion**    | T(x)       | τ(x)       | θ(x)         |



lb.behaviour(axis, points).quantity.modifier

lb.shear(y, 100).force.max.at(2.5)

behaviours are:

- shear
- bending
- axial
- torsion

quantities are:

- stress
- displacement

modifiers will be:
- max
- mean
- min
- range
- at(x)

not using a modifier will return the entire array
all actions will have the same modifier options

axis does not need to be an input for bending or torsion
e.g.
lb.shear(y, 100).stress.range

also the loaded beam will have a deflection method (that will take axis and points as inputs)
e.g.
lb.deflection(z, 100).max

## von Mises

i want to have von mises as follows

lb.von_mises(points).modifier

