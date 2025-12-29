from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt

from beamy.results.member_result import MemberResult


def plot_member_axial(member_result: MemberResult, save_path: Optional[str] = None, show: bool = False) -> None:
    fig, ax = plt.subplots()
    stations = [0.0, 1.0]
    axial = [member_result.axial_force, member_result.axial_force]
    ax.plot(stations, axial, label="Axial (N)", color="darkred")
    ax.fill_between(stations, axial, color="salmon", alpha=0.3)
    ax.set_xlabel("Station (0-1)")
    ax.set_ylabel("Axial force")
    ax.set_title(f"Member {member_result.member_id} axial")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend()
    plt.tight_layout()

    if save_path:
        if not save_path.lower().endswith(".svg"):
            save_path = f"{save_path}.svg"
        plt.savefig(save_path, format="svg")
    if show and not save_path:
        plt.show()
    plt.close(fig)

