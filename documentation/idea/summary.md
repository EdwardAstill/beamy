### Beamy (project summary)

Beamy is a lightweight structural analysis package for building simple frame/beam-style models, applying loads and supports, and getting back the key results an engineer cares about—deflections, reactions, and internal member forces. The goal is to keep the workflow clear and predictable: you describe the structure, describe what’s applied to it, then run an analysis to see how it behaves.

### What you use it for

- **Define a structure**: create nodes (points in space) and connect them with members (beams/trusses/cables).
- **Apply actions**: specify loads (forces/moments) and supports (restraints) without mixing them into the core geometry (supports are per-load-case boundary conditions).
- **Run analysis**: solve for displacements and reactions, and recover forces along members.
- **Review results**: inspect member diagrams and overall structural response, and optionally visualize the model and results.

### Design goals

- **Engineer-friendly inputs**: models map closely to how engineers describe structures.
- **Separation of concerns**: geometry/topology, loading, analysis settings, and results are kept distinct.
- **Extensible foundation**: structured so more capabilities (design checks, additional solvers, richer connections) can be added over time.