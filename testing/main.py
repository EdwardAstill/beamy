from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from sectiony.library import rhs

# Add the src directory to the path so we can import beamy as a package
SRC_PATH = Path(__file__).resolve().parents[1] / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from beamy import LoadCase, LoadedMember, Material, MemberPointForce, plot_beam_diagram


def main() -> None:
    mat = Material(name="Steel", E=200e9, G=80e9)
    sec = rhs(b=0.1, h=0.2, t=0.005, r=0.0)

    L = 1.0  # m

    lc = LoadCase(name="Case 1")
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=L / 2.0,
            force=np.array([0.0, 0.0, -10_000.0]),
            coords="global",
            position_type="absolute",
        )
    )

    lm = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=sec,
        material=mat,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start="111100",
        support_end="011000",
        load_case=lc,
    )

    profile = lm.member_demand().actions(points=401)
    print("Max |Vz| =", profile.shear_z.abs_max)
    print("Max |My| =", profile.bending_y.abs_max)

    # Plot (interactive)
    plot_beam_diagram(lm, plot_stress=False, plot_section=True, save_path=None)


if __name__ == "__main__":
    main()

