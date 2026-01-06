# Analysis Module Reference

The `analysis/` module handles all FE mechanics, meshing, assembly, solving, and recovery. It is the core "computational engine" of Beamy.

## `analysis/model.py`

Internal solver model structures created inside `analyze()`. These are "dumb containers" separate from the user-facing `model/` classes.

*   `AnalysisNode(id, xyz)`: Internal node representation.
*   `Element(id, type, node_i, node_j, ...)`: Represents a mesh element (segment of a member).
*   **Mappings**:
    *   `member_to_elements`: Map `MemberId` to list of `ElementId`.
    *   `station_node_map`: Map `(MemberId, s)` to `AnalysisNodeId`.
    *   `element_state`: Stores history for second-order analysis (axial force, corotational frame).

## `analysis/settings.py`

Defines the configuration for an analysis run.

*   `AnalysisSettings`:
    *   `dimension`: `"3d"` (future `"2d"`).
    *   `solver`: `"linear"` | `"second_order"`.
    *   `formulation`: `"linear"` | `"pdelta"` | `"corotational"`.
    *   **Tolerances**: Max iterations, load steps, convergence criteria.
    *   **Meshing**: Options for automatic station insertion.

## `analysis/mesh.py`

Converts `Frame + LoadCase` into an `AnalysisModel`.

*   **Automated Station Discovery (Station Manager)**:
    1.  **Collect** all station ($s$) values from:
        *   `Member.stations` (user-defined).
        *   `PointLoad` locations.
        *   `Support` locations.
        *   `Connection` (`EndToStation`) targets.
    2.  **Deduplicate** using `math.isclose` with tolerance (e.g., `merge_tol/L`).
    3.  **Sort** stations: `[0.0, ..., 1.0]`.
    4.  **P-Delta Refinement**: If `formulation="pdelta"`, forces insertion of at least one mid-point station ($s=0.5$) to capture member curvature effects.
*   **Process**:
    *   Create `AnalysisNode`s at all stations.
    *   Split members into `Element`s between nodes.
    *   Apply end-to-station connectivity by wiring endpoints to the inserted station nodes.
    *   **Preserve Mapping**: Store `station_node_map` in `FrameResult` to allow result queries at specific physical locations (e.g., "moment at support").

## `analysis/transformations.py`

Handles local/global coordinate transformations.

*   **Functions**: Direction cosines, rotation matrices, element transformation `T`.
*   **Vertical Member Handling**:
    *   Detects if a member is parallel to the global up vector.
    *   Falls back to "Global North" or "Global East" reference vectors to resolve orientation ambiguity.

## `analysis/elements/`

Element implementations and factory.

*   `make_element(element, frame, settings)`: Factory to create element instances.
*   `beam3d.py`: Core 3D beam element.
    *   `k_material(...)`: Local stiffness matrix.
    *   `k_geometric(...)`: Geometric stiffness (for P-Delta).
    *   **Static Condensation**:
        *   Handles member end releases at the **Element Level**.
        *   Partitions $K_e$ and condenses out released DOFs.
        *   Stores **Condensation Matrix** ($K_{cc}^{-1} K_{cr}$) to allow `recovery.py` to calculate internal displacements at released ends ($u_r$).
    *   **Load Conversion**:
        *   Calculates **Fixed-End Actions** for member loads based on element shape functions.
*   `truss3d.py`: Axial-only stiffness.
*   `cable3d.py`: Tension-only behavior (nonlinear iteration logic).

## `analysis/assembly.py`

Builds the global system of equations.

*   **DOF Numbering**: Maps `AnalysisNode` DOFs to global equation indices.
*   **Assembly**: Aggregates global `K` and `F`.
    *   Asks elements for **Fixed-End Actions** when building the force vector $F$.
*   **Releases**: Handled by elements returning condensed matrices.
*   **Supports**: Applies restraints, prescribed displacements, and springs.

## `analysis/solver_linear.py`

Solves linear systems.

*   `u = solve(K, F)`
*   Computes reactions.
*   Strictly mathematical; no plotting or result formatting.

## `analysis/solver_second_order.py`

Driver for second-order/nonlinear analysis.

*   **Process**:
    *   Load stepping.
    *   Newton/Modified Newton iterations.
    *   Assembles tangent stiffness $K_t = K_m + K_g$.
    *   Updates element state (axial forces) for $K_g$.
    *   Checks convergence.

## `analysis/recovery.py`

Post-processes solution `u` into engineering outputs.

*   Converts `u` to element end forces.
*   Computes internal force diagrams at requested stations.
*   Computes stresses (axial + bending).
*   Aggregates element results back to parent `MemberResult`.

