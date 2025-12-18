"""Smoke test for locally-defined loads in frame analysis.

This is intentionally minimal and prints key values for sanity.
"""

from __future__ import annotations

import numpy as np

from beamy.frame import FrameBuilder
from beamy.core.material import Material
from sectiony.library import rhs
from beamy.core.loads import FrameLoadCase
from beamy.frame.analysis import LoadedFrame


def main() -> None:
    # Simple cantilever member along global X.
    material = Material(name="Steel", E=200e9, G=79.3e9)
    section = rhs(b=0.1, h=0.2, t=0.005, r=0.0)

    b = FrameBuilder()
    b.add("M1", (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), section=section, material=material, orientation=(0, 0, 1))
    # Fix start node (builder will create N0 at x=0 and N1 at x=2)
    b.support_at((0.0, 0.0, 0.0), "111111")
    frame = b.build()

    lc = FrameLoadCase("local-loads")

    # Apply a local member point moment about local z at midspan.
    # local z for this member (x along +X, orientation +Z => local y ~ +Z, local z ~ -Y)
    lc.add_member_point_moment("M1", position=1.0, moment=np.array([0.0, 0.0, 1_000.0]), coords="local")

    # Apply a nodal force at the free end, defined in member local coordinates.
    # Use reference_member_id to define the local axes.
    free_node = frame.get_member("M1").end_node_id
    lc.add_nodal_force(free_node, np.array([0.0, 100.0, 0.0]), coords="local", reference_member_id="M1")

    loaded = LoadedFrame(frame, lc)

    # Sanity prints
    r0 = loaded.reactions[frame.get_member("M1").start_node_id]
    print("Reactions at fixed end [Fx,Fy,Fz,Mx,My,Mz] =", np.round(r0, 6))

    beams = loaded.to_loaded_beams()
    lb = beams["M1"]
    # Check that beam loadcase has the point moment we injected
    print("Beam moments (x, [T,My,Mz]) =", [(m.x, np.round(m.moment, 6)) for m in lb.loads.moments])


if __name__ == "__main__":
    main()
