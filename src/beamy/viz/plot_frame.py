from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from beamy.model.frame import Frame


def plot_frame_undeformed(frame: Frame, save_path: Optional[str] = None, show: bool = False) -> None:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # plot members
    for member in frame.members.values():
        node_i = frame.nodes[member.end_i]
        node_j = frame.nodes[member.end_j]
        xs = [node_i.xyz[0], node_j.xyz[0]]
        ys = [node_i.xyz[1], node_j.xyz[1]]
        zs = [node_i.xyz[2], node_j.xyz[2]]
        ax.plot(xs, ys, zs, color="steelblue", linewidth=2)

    # plot nodes
    for node in frame.nodes.values():
        ax.scatter(node.xyz[0], node.xyz[1], node.xyz[2], color="black", s=12)
        ax.text(node.xyz[0], node.xyz[1], node.xyz[2], f"{node.id}", fontsize=8)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=20, azim=45)
    plt.tight_layout()

    if save_path:
        if not save_path.lower().endswith(".svg"):
            save_path = f"{save_path}.svg"
        plt.savefig(save_path, format="svg")
    if show and not save_path:
        plt.show()
    plt.close(fig)

