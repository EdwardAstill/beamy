# main.py
import sys
from pathlib import Path

# Add the src directory to the path so we can import beamy as a package
src_path = Path(__file__).parent.parent
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import numpy as np
from beamy.beam import Beam1D, Material, Section
from beamy.loads import LoadCase, PointLoad
from beamy.analysis import analyze_beam_simple_point_load

def main():
    mat = Material(name="Steel", E=200e9, G=80e9)
    sec = Section(name="Test", A=0.01, Iy=1e-6, Iz=1e-6, J=1e-6, y_max=0.1, z_max=0.1)
    L = 1.0  # 1 m beam

    beam = Beam1D(
        L=L,
        material=mat,
        section=sec,
        support_left="111000",  # Pinned: Ux, Uy, Uz constrained; Rx, Ry, Rz free
        support_right="111000",  # Pinned: Ux, Uy, Uz constrained; Rx, Ry, Rz free
    )

    # 10 kN downward at midspan
    P = -10_000.0  # negative in local z
    load = PointLoad(
        x=0.5,
        force=np.array([0.0, 0.0, P]),
    )

    lc = LoadCase(name="Case 1")
    lc.add_point_load(load)

    result = analyze_beam_simple_point_load(beam, lc)

    print("Reactions:", result.reactions)
    for x in np.linspace(0.0, L, 11):
        print(f"x={x:.2f} m: Vz={result.Vz(x):.2f}, Mz={result.Mz(x):.2f}, w={result.w(x):.6e}")

if __name__ == "__main__":
    main()
