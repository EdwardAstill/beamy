# Beamy Examples

This directory contains example scripts demonstrating various features of Beamy.

## Examples

### `simple_beam.py`
A basic simply supported beam with a point load at mid-span. Demonstrates:
- Basic beam setup
- Point loads
- Simple visualization

### `cantilever.py`
A cantilever beam with a point load at the free end. Demonstrates:
- Fixed support conditions
- End loading
- Stress visualization

### `ibeam.py`
An I-beam with multiple supports and various load types. Demonstrates:
- Multiple supports
- Point forces and distributed forces
- Finding critical stress locations
- Section stress plotting

### `multi_load.py`
A beam with multiple load types including:
- Point forces
- Distributed forces
- Moments
- Multiple supports

### `stress_analysis.py`
Comprehensive stress analysis example showing:
- Finding maximum stress locations
- Plotting section stresses at critical locations
- Comparing different stress types (bending, shear, von Mises)
- Multiple visualization outputs

### `different_sections.py`
Comparison of different cross-section types:
- I-beam
- Rectangular Hollow Section (RHS)
- Shows how section shape affects beam behavior

### `line_plots.py`
Demonstrates how to generate 2D line plots for:
- Shear Force Diagrams
- Bending Moment Diagrams
- Deflection Diagrams
- Axial/Torsion Diagrams

### `support_reactions.py`
Demonstrates how to access numerical support reaction values:
- Continuous beam with multiple supports
- Complex multi-axial loading
- Accessing reaction forces (Fx, Fy, Fz) and moments (Mx, My, Mz) programmatically

### `support_plot.py`
Demonstrates how to visualize the beam supports in 2D.

## Running Examples

All examples save their plots to the `gallery/` directory. 

**Note:** Make sure Beamy is installed in your Python environment. If running from the project root, you may need to install it in development mode:

```bash
pip install -e .
```

To run an example:

```bash
# From the project root
python examples/simple_beam.py

# Or from the examples directory
cd examples
python simple_beam.py
```

The plots will be saved as PNG files in the `gallery/` directory with descriptive names.

## Gallery

After running the examples, check the `gallery/` directory to see the generated visualizations.

