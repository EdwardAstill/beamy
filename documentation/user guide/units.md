# Units

Beamy is **unit-agnostic**. You can use any consistent system of units (SI, US customary, etc.), as long as all quantities use the same base units throughout your analysis.

## Unit Consistency

**Critical:** All inputs for a single analysis must use the same unit system. Mixing units will produce incorrect results.

## Unit Relationships

The following relationships must be maintained:

- **Stress** = Force / Length²
- **Moment** = Force × Length
- **Young's Modulus** = Force / Length² (same units as stress)
- **Shear Modulus** = Force / Length² (same units as stress)
- **Area** = Length²
- **Second Moment of Area** = Length⁴
- **Distributed Load** = Force / Length

## Example Unit Systems

### SI Units
- Length: meters (m)
- Force: Newtons (N)
- Stress/Modulus: Pascals (Pa = N/m²)
- Moment: Newton-meters (N·m)

### US Customary
- Length: inches (in) or feet (ft)
- Force: pounds-force (lbf) or kips (kip = 1000 lbf)
- Stress/Modulus: psi (lbf/in²) or ksi (kip/in²)
- Moment: lbf·in, lbf·ft, or kip·in, kip·ft

### Metric (Alternative)
- Length: millimeters (mm)
- Force: Newtons (N)
- Stress/Modulus: Megapascals (MPa = N/mm²)
- Moment: Newton-millimeters (N·mm)

## Best Practice

Choose one unit system and use it consistently for all inputs:
- All lengths in the same unit
- All forces in the same unit
- All material properties (E, G) in force/length²
- All section properties in length² and length⁴

The output results will automatically use the same units as your inputs.