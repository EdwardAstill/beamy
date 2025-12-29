"""
Different Section Types Example (frame-first API)

Demonstrates how different cross-sections change response under the same load case.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from sectiony.library import i_section, rhs

from beamy import LoadCase, LoadedMember, Material, MemberPointForce, plot_beam_diagram


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    gallery_dir = project_root / "gallery" / "section"
    gallery_dir.mkdir(parents=True, exist_ok=True)

    steel = Material(name="Steel", E=200e9, G=80e9)
    L = 4.0

    # Common loading: 10 kN downward at midspan (global -Z)
    load_case = LoadCase(name="Point Load")
    load_case.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=L / 2.0,
            force=np.array([0.0, 0.0, -10_000.0]),
            coords="global",
            position_type="absolute",
        )
    )

    # Common supports: constrain translations + Rx at start; constrain Uy/Uz at end
    support_start = "111100"
    support_end = "011000"

    # 1) I-Beam
    print("=== I-Beam ===")
    section_ibeam = i_section(d=0.2, b=0.1, tf=0.01, tw=0.006, r=0.008)
    lb_ibeam = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section_ibeam,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start=support_start,
        support_end=support_end,
        load_case=load_case,
    )

    plot_beam_diagram(
        lb_ibeam,
        plot_stress=True,
        plot_section=True,
        save_path=str(gallery_dir / "section_ibeam.svg"),
    )

    # 2) RHS
    print("\n=== Rectangular Hollow Section ===")
    section_rect = rhs(b=0.1, h=0.2, t=0.005, r=0.0)
    lb_rect = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section_rect,
        material=steel,
        orientation=np.array([0.0, 1.0, 0.0]),
        support_start=support_start,
        support_end=support_end,
        load_case=load_case,
    )

    plot_beam_diagram(
        lb_rect,
        plot_stress=True,
        plot_section=True,
        save_path=str(gallery_dir / "section_rhs.svg"),
    )

    print("\nPlots saved:")
    print(f"  - {gallery_dir / 'section_ibeam.svg'}")
    print(f"  - {gallery_dir / 'section_rhs.svg'}")


if __name__ == "__main__":
    main()