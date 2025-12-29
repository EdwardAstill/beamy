from __future__ import annotations

import sys
from pathlib import Path

from sectiony.library import rhs

# Add src to path for development
SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from beamy import plot_section


def test_section_plot() -> None:
    sec = rhs(b=0.1, h=0.2, t=0.005, r=0.0)
    plot_section(sec, show=True)


if __name__ == "__main__":
    test_section_plot()

