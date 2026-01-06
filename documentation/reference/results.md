# Results Module Reference

The `results/` module provides structured outputs for analysis runs. Result objects are separated from model objects to allow comparison of multiple load cases.

## `results/frame_result.py`

*   `FrameResult`: A standalone snapshot of the entire structure's response.
    *   **Indexed by Element IDs**: Keys correspond to elements at the time of the solve.
    *   `u`: Nodal displacements.
    *   `reactions`: Support reactions.
    *   `element_end_forces`: Forces at element ends.
    *   `member_diagrams`: Aggregated internal force diagrams for members.
    *   `station_node_map`: Mapping of `(MemberId, s)` to internal node IDs for location-based querying.
    *   **Metadata**: Convergence status, step history, iteration counts.

## `results/member_result.py`

*   `MemberResult`: Detailed response for a single member.
    *   **Diagrams**: Axial (N), Shear (V), Moment (M), Torsion (T).
    *   **Deflections**: Deformed shape along the span.
    *   **Stresses**: Stress arrays at stations.
    *   **Mapping**: Access values by station `s`.
    *   Typically retrieved via `FrameResult.member(member_id)`.

