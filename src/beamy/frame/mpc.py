from __future__ import annotations
from dataclasses import dataclass
from typing import Literal


@dataclass
class MPC:
    """
    Multi-point constraint tying a free node to an interior point along a member.

    The constraint enforces (selected) DOFs of the node to match the interpolated
    DOFs of the member at a relative position xi (0 = start, 1 = end).

    dofs is a 6-char string "111111" selecting [UX, UY, UZ, RX, RY, RZ].
    """
    member_id: str
    xi: float  # 0..1 along the member
    node_id: str
    dofs: Literal["000000", "111111"] | str = "111111"

    def __post_init__(self) -> None:
        if not (0.0 <= self.xi <= 1.0):
            raise ValueError(f"MPC xi must be in [0,1], got {self.xi}")
        if not isinstance(self.dofs, str) or len(self.dofs) != 6 or not all(c in "01" for c in self.dofs):
            raise ValueError(f"MPC dofs must be 6-char 0/1 string, got '{self.dofs}'")
