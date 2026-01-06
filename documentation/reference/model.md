# Model Module Reference

The `model/` module contains the user-facing structural object definitions. These classes are primarily data containers with validation and convenience wrappers.

## `model/node.py`

*   `Node(id, xyz, meta)`: Represents a spatial point.
    *   Basic validation of coordinates.
    *   No merging logic (handled by `Frame`).

## `model/member.py`

*   `Member`: Represents a physical structural member.
    *   **Topology**: `end_i`, `end_j` (Node IDs).
    *   **Type**: `kind` ("beam", "truss", "cable").
    *   **Properties**: `material_id`, `section_id` (Refs to Frame registries).
    *   **Orientation**: `local_z`, `roll_angle`, or `ref_node`.
    *   **Releases**: `releases_i`, `releases_j` (`ReleaseMask6`).
    *   **Stations**: `stations` set (mesh seeds, 0..1).
    *   **Design Properties**: `Lx`, `Ly` (unbraced), `Kx`, `Ky`, `Kz`, `Cb`.
        *   *Validation*: Missing values (e.g., $L_b$) should be checked during `Frame.validate()`.

## `model/frame.py`

*   `Frame`: The root container for the structural model.
    *   **Storage**: `nodes`, `members`, `connections`.
    *   **Registries**:
        *   `materials`: `dict[str, Material]`
        *   `sections`: `SectionRegistry` (Smart lookup for standard shapes like "W14x90").
    *   **Methods**:
        *   `add_node()`, `add_member()`.
        *   `connect_end_to_member()`: Creates `EndToStation` connection.
        *   `build()`: Validates model, materializes split requirements from connections/loads.
        *   `analyze(loadcase, settings)`: Delegates to `analysis/` via local import.
        *   `validate()`: Checks IDs, properties, and station ranges.
            *   Ensures all `member.section_id` references exist in `SectionRegistry`.
            *   Validates member design properties (e.g., defaults for unbraced length).

## `model/dof.py`

*   **DOF Order (3D)**: `[UX, UY, UZ, RX, RY, RZ]`.
*   **Mask Types**:
    *   `ReleaseMask6`: Member end releases (True = Released/Free).
    *   `RestraintMask6`: Supports (True = Restrained/Fixed).

## `model/connections.py`

Typed connection definitions used to describe topology.

*   `EndRef(member_id, end)`: References a specific end of a member.
*   `StationRef(member_id, s)`: References a point along a member.
*   `EndToStation(end, target, dofs)`: Connects a member end to a point on another member.

