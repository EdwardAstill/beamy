import sys
from pathlib import Path
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import numpy as np
from sectiony.library import rhs

from beamy import LoadCase, LoadedMember, Material, MemberPointForce, plot_beam_diagram

if __name__ == "__main__":
    mat = Material(name="Test", E=200e9, G=80e9)
    sec = rhs(b=0.1, h=0.2, t=0.005, r=0.0)
    L = 5.0
    
    # Create a test load case
    lc = LoadCase(name="Test Case")
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=2.5,
            force=np.array([0.0, -1000.0, 500.0]),
            coords="global",
            position_type="absolute",
        )
    )
    lc.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=5.0,
            force=np.array([-100.0, 0.0, 0.0]),
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
        support_start="111111",
        support_end="111111",
        load_case=lc,
    )

    plot_beam_diagram(lm, plot_stress=False, plot_section=True, save_path=None)
