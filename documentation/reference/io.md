# I/O Module Reference

The `io/` module handles import and export of model data and results.

## `io/json.py`

JSON serialization and deserialization.

*   **Supported Types**:
    *   `Frame`
    *   `LoadCase`
    *   `AnalysisSettings`
*   **Features**:
    *   Stable schema with versioning.
    *   Handles dataclasses and custom types (e.g., enums).

## `io/csv.py`

Export of analysis results to tabular format (CSV).

*   **Exports**:
    *   Nodal displacements.
    *   Support reactions.
    *   Member force and stress summaries.
*   **Scope**: Data export only; no plotting or solving logic.

