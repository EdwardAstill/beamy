# Beamy

Beamy is a lightweight Python package for 1D beam analysis, capable of handling static loads using Euler-Bernoulli beam theory. It supports axial, torsional, and transverse (shear/bending) analysis with 6 degrees of freedom per node.

## Features

- **Unit-agnostic**: Use any consistent unit system (e.g., SI or Imperial).
- **Beam Theory**: Euler-Bernoulli beam formulation.
- **Loading**:
    - Point forces and moments (eccentric loads supported).
    - Distributed loads (linearly varying).
    - Combined axial, torsional, and bending loads.
- **Supports**: Fully customizable 6-DOF support conditions (fixed, pinned, roller, etc.).
- **Sections**: Integration with `sectiony` for complex cross-section properties.
- **Analysis**:
    - Reaction forces and moments.
    - Internal force/moment diagrams (Shear, Bending, Axial, Torsion).
    - Deflection and rotation profiles.
    - Von Mises stress distribution along the beam.
- **Visualization**:
    - 3D interactive beam diagrams with loads and supports.
    - Stress distribution plots on cross-sections.

## Installation

```bash
pip install .
```

*Note: Beamy requires `numpy`, `matplotlib`, and `sectiony`.*

## Quick Start

Here is a simple example of a simply supported beam with a point load.

```python
import numpy as np
from beamy import Beam1D, Material, Support, LoadCase, PointForce, LoadedBeam
from sectiony.library import i_section

# 1. Define Material and Section
steel = Material(name="Steel", E=200e9, G=80e9)
section = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.0)

# 2. Create Beam
# 5m long beam, Pinned at x=0, Roller at x=5
L = 5.0
supports = [
    Support(x=0.0, type="111100"), # Pinned (x,y,z translation fixed, rotation free)
    Support(x=L, type="011000")    # Roller (y,z translation fixed, x free, rotation free)
]
beam = Beam1D(L=L, material=steel, section=section, supports=supports)

# 3. Apply Loads
# 10kN downward force at mid-span
load = PointForce(
    point=np.array([2.5, 0.0, 0.0]),
    force=np.array([0.0, 0.0, -10000.0])
)
case = LoadCase(name="Design Load")
case.add_point_force(load)

# 4. Solve
lb = LoadedBeam(beam, case)

# 5. Get Results
print(f"Max Deflection: {lb.deflection('z').abs_max:.6f} m")
print(f"Max Bending Moment: {lb.bending('z').action.abs_max:.2f} Nm")

# 6. Visualize
# Plot 3D beam diagram with stress coloring
lb.plot(plot_stress=True)
```

## Documentation

For more detailed usage and examples, please refer to the [Documentation](documentation/):

- [User Guide](documentation/user%20guide/user_guide.md): Comprehensive guide on setting up beams, loads, and interpreting results.
- [Reference](documentation/reference/): detailed API reference.

## License

MIT

