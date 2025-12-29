from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Union


@dataclass
class MemberResult:
    member_id: Union[int, str]
    kind: str
    axial_force: float
    end_forces_local: Optional[Tuple[float, ...]] = None

    def utilization(self) -> float:
        # placeholder until design checks exist
        return 0.0

