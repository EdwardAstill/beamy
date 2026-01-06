# Loads Module Reference

The `loads/` module defines inputs for analysis cases. These classes define *applied* loads and boundary conditions; they contain no solver math.

## `loads/loadcase.py`

*   `LoadCase`: Container for all loads and supports in a specific scenario.
    *   `name`: str
    *   `nodal_loads`: list[`NodalLoad`]
    *   `member_loads`: list[`MemberLoad`]
    *   `supports`: list[`Support`]
    *   `prescribed_displacements`: list
    *   `springs`: list

## `loads/nodal_loads.py`

*   `NodalLoad`: Force/moment applied to a node.
    *   `node_id`: ID of target node.
    *   `forces`: $(F_x, F_y, F_z)$
    *   `moments`: $(M_x, M_y, M_z)$
    *   `coord_sys`: `"global"` | `"local"`

## `loads/member_loads.py`

Defines loads applied along member spans.
*   **Note**: The logic to convert these to nodal equivalents (Fixed-End Actions) resides in `analysis/elements/beam3d.py` as it depends on element shape functions.

*   `MemberPointLoad`: Concentrated force at station $s$.
    *   `member_id`, `s`, `vector`, `coord_sys`.
*   `MemberPointMoment`: Concentrated moment at station $s$.
*   `MemberDistributedLoad`: Line load from $s_0$ to $s_1$.
    *   `w0`, `w1`: Intensities at start/end.

## `loads/supports.py`

Defines boundary conditions and constraints.

*   `Support`: Restrains specific degrees of freedom.
    *   `target`: `NodeId` or `StationRef`.
    *   `restrained_dofs`: `RestraintMask6` (True = restrained/fixed).
*   `PrescribedDisplacement`: Enforces non-zero displacement at a DOF.
*   `Spring`: Elastic support.
    *   `k_trans`, `k_rot`.

