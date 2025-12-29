from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Tuple

from beamy.results.frame_result import FrameResult

Disp6 = Tuple[float, float, float, float, float, float]


def export_node_displacements(result: FrameResult, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "ux", "uy", "uz", "rx", "ry", "rz"])
        for node_id, disp in result.node_displacements.items():
            writer.writerow([node_id, disp[0], disp[1], disp[2], disp[3], disp[4], disp[5]])


def export_reactions(result: FrameResult, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "fx", "fy", "fz", "mx", "my", "mz"])
        for node_id, reac in result.reactions.items():
            writer.writerow([node_id, reac[0], reac[1], reac[2], reac[3], reac[4], reac[5]])


def export_member_axial_forces(result: FrameResult, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["member_id", "kind", "axial_force"])
        for member_id, mr in result.member_results.items():
            writer.writerow([member_id, mr.kind, mr.axial_force])

