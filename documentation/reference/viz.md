# Viz Module Reference

The `viz/` module provides plotting capabilities using `matplotlib`. It is read-only and contains no solver logic.

## `viz/plot_frame.py`

Functions for plotting the global structure.

*   **Inputs**: `Frame` object, optional style configurations.
*   **Features**:
    *   Plots nodes and members.
    *   Labels for IDs.
    *   Markers for supports and releases.
    *   visualization of local axes.

## `viz/plot_member.py`

Functions for detailed member-level plotting.

*   **Inputs**: `Frame`, `MemberResult` or `FrameResult`.
*   **Features**:
    *   Internal force diagrams (N, V, M, T).
    *   Stress distribution plots.
    *   Deflection diagrams.

