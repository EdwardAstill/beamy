from __future__ import annotations

from pathlib import Path

import numpy as np
from sectiony.library import i as i_section

from beamy import (
    LoadCase,
    LoadedMember,
    Material,
    MemberDistributedForce,
    MemberPointForce,
    MemberPointMoment,
    plot_analysis_results,
)


def main() -> None:
    # 1. Setup member (units: mm, N)
    project_root = Path(__file__).resolve().parents[2]
    gallery_dir = project_root / "gallery" / "plots"
    gallery_dir.mkdir(parents=True, exist_ok=True)

    section = i_section(d=200, b=100, tf=10, tw=6, r=8)
    steel = Material(name="Steel", E=210e9, G=80e9)

    L = 3000.0  # mm

    # 2. Setup loads
    load_case = LoadCase(name="Complex Loading")

    # Vertical point force at x=1000 mm (global +Y is "up" here; use negative for down)
    load_case.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=1000.0,
            force=np.array([0.0, -5000.0, 0.0]),
            coords="global",
            position_type="absolute",
        )
    )

    # Horizontal point force at x=2000 mm (global +Z)
    load_case.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=2000.0,
            force=np.array([0.0, 0.0, 2000.0]),
            coords="global",
            position_type="absolute",
        )
    )

    # Distributed force from x=1500 to 2500 mm (global -Y)
    load_case.member_distributed_forces.append(
        MemberDistributedForce(
            member_id="M1",
            start_position=1500.0,
            end_position=2500.0,
            start_force=np.array([0.0, -2000.0, 0.0]),
            end_force=np.array([0.0, -2000.0, 0.0]),
            coords="global",
        )
    )

    # Axial force at x=L (compression in global -X)
    load_case.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=L,
            force=np.array([-10000.0, 0.0, 0.0]),
            coords="global",
            position_type="absolute",
        )
    )

    # Torsion about global X at x=1500 mm
    load_case.member_point_moments.append(
        MemberPointMoment(
            member_id="M1",
            position=1500.0,
            moment=np.array([500.0, 0.0, 0.0]),
            coords="global",
            position_type="absolute",
        )
    )

    # 3. Analyze (supports: pinned-like at start, roller-like at end)
    lm = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 0.0, 1.0]),
        support_start="111100",
        support_end="011100",
        load_case=load_case,
    )

    # 4. Plot results (saved as .svg)
    out_path = gallery_dir / "line_plots_example.svg"
    plot_analysis_results(
        lm,
        save_path=str(out_path),
        show=False,
        points=201,
        units={"length": "mm", "force": "N", "moment": "N·mm", "deflection": "mm"},
    )
    print(f"Plots saved to {out_path}")


if __name__ == "__main__":
    main()