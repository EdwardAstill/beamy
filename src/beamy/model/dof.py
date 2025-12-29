from __future__ import annotations

from typing import Iterable, Literal, Tuple

DOF_ORDER: Tuple[str, ...] = ("UX", "UY", "UZ", "RX", "RY", "RZ")

DofLabel = Literal["UX", "UY", "UZ", "RX", "RY", "RZ"]
ReleaseMask6 = Tuple[bool, bool, bool, bool, bool, bool]
RestraintMask6 = Tuple[bool, bool, bool, bool, bool, bool]


def dof_index(label: DofLabel) -> int:
    index = DOF_ORDER.index(label)
    return int(index)


def make_release_mask(mask: Iterable[bool]) -> ReleaseMask6:
    values = tuple(bool(value) for value in mask)
    _ensure_len(values, "release mask")
    return values  # type: ignore[return-value]


def make_restraint_mask(mask: Iterable[bool]) -> RestraintMask6:
    values = tuple(bool(value) for value in mask)
    _ensure_len(values, "restraint mask")
    return values  # type: ignore[return-value]


def _ensure_len(values: Tuple[bool, ...], label: str) -> None:
    if len(values) != 6:
        raise ValueError(f"{label} must contain 6 entries")

