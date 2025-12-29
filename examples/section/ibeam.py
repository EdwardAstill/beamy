from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sectiony.library import i as i_section

from beamy import (
    LoadCase,
    LoadedMember,
    Material,
    MemberDistributedForce,
    MemberPointForce,
    MemberPointSupport,
    StressPlotter,
    plot_beam_diagram,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    # Units: mm, N
    section = i_section(d=200, b=100, tf=10, tw=6, r=8)
    steel = Material(name="Steel", E=210e9, G=80e9)

    L = 3000.0  # mm

    # Loading
    loads = LoadCase(name="Test Load")
    loads.member_point_forces.append(
        MemberPointForce(
            member_id="M1",
            position=1000.0,
            force=np.array([0.0, -5000.0, 0.0]),
            coords="global",
            position_type="absolute",
        )
    )
    loads.member_distributed_forces.append(
        MemberDistributedForce(
            member_id="M1",
            start_position=1300.0,
            end_position=2500.0,
            start_force=np.array([0.0, -5.0, 0.0]),
            end_force=np.array([0.0, -1.0, 0.0]),
            coords="global",
        )
    )

    # Supports: start/end + an intermediate support at x=1300 mm
    lb = LoadedMember(
        id="M1",
        start=np.array([0.0, 0.0, 0.0]),
        end=np.array([L, 0.0, 0.0]),
        section=section,
        material=steel,
        orientation=np.array([0.0, 0.0, 1.0]),
        support_start="111100",
        support_end="011100",
        point_supports=[MemberPointSupport(position=1300.0, support="111111")],
        load_case=loads,
    )

    # Find the location of maximum Von Mises stress
    vm_results = lb.von_mises(points=401)
    max_vm_idx = int(np.argmax(vm_results._values))
    max_vm_x = float(vm_results._x[max_vm_idx])
    max_vm_val = float(vm_results._values[max_vm_idx])

    print(f"Max Von Mises Stress: {max_vm_val/1e6:.2f} MPa at x = {max_vm_x:.2f} mm")

    gallery_dir = PROJECT_ROOT / "gallery" / "section"
    gallery_dir.mkdir(parents=True, exist_ok=True)

    # Save beam diagram
    plot_beam_diagram(
        lb,
        plot_stress=True,
        plot_section=True,
        save_path=str(gallery_dir / "ibeam_3d_diagram.svg"),
    )

    # Plot section stress at the critical location
    sp = StressPlotter(lb)
    sp.plot_stress_at(
        x_pos=max_vm_x,
        stress_type="von_mises",
        title=f"Max Von Mises Stress Section (x={max_vm_x:.0f}mm)",
        cmap="plasma",
        show=False,
    )
    plt.savefig(gallery_dir / "ibeam_section_stress.svg", bbox_inches="tight", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()