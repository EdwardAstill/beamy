from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

if TYPE_CHECKING:
    from .section import Section

def plot_section(section: Section, ax: Optional[plt.Axes] = None, show: bool = True) -> Optional[plt.Axes]:
    """
    Plot the cross-section geometry.
    
    Args:
        section: The section to plot
        ax: Optional matplotlib axes to plot on. If None, creates a new figure.
        show: Whether to call plt.show() at the end.
    
    Returns:
        The axes object used for plotting.
    """
    if not section.geometry or not section.geometry.shapes:
        print("No geometry defined for this section.")
        return None

    if ax is None:
        fig, ax = plt.subplots()
    
    # Shapes are polygons of (y, z) points
    # We want to plot as (z, y) so z is horizontal, y is vertical
    
    solids = [s for s in section.geometry.shapes if not s.hollow]
    hollows = [s for s in section.geometry.shapes if s.hollow]
    
    all_z = []
    all_y = []
    
    # Plot solids
    for s in solids:
        # Convert (y, z) to (z, y)
        vertices = [(p[1], p[0]) for p in s.points]
        zs = [p[1] for p in s.points]
        ys = [p[0] for p in s.points]
        all_z.extend(zs)
        all_y.extend(ys)
        
        poly = Polygon(vertices, closed=True, facecolor='silver', edgecolor='black', alpha=0.8, label='Solid')
        ax.add_patch(poly)
        
    # Plot hollows
    for h in hollows:
        vertices = [(p[1], p[0]) for p in h.points]
        zs = [p[1] for p in h.points]
        ys = [p[0] for p in h.points]
        all_z.extend(zs)
        all_y.extend(ys)
        
        poly = Polygon(vertices, closed=True, facecolor='white', edgecolor='black', linestyle='--', alpha=1.0, label='Hollow')
        ax.add_patch(poly)
        
    # Set limits and aspect
    if all_z and all_y:
        z_min, z_max = min(all_z), max(all_z)
        y_min, y_max = min(all_y), max(all_y)
        
        dz = z_max - z_min
        dy = y_max - y_min
        
        # Handle case where section is a line or point
        if dz == 0: dz = 1.0
        if dy == 0: dy = 1.0
        
        padding_z = dz * 0.1
        padding_y = dy * 0.1
        
        ax.set_xlim(z_min - padding_z, z_max + padding_z)
        ax.set_ylim(y_min - padding_y, y_max + padding_y)
        ax.set_aspect('equal')
    
    ax.set_xlabel('z (Horizontal)')
    ax.set_ylabel('y (Vertical)')
    ax.set_title(f"Section: {section.name}")
    ax.grid(True, linestyle=':', alpha=0.6)

    if show:
        plt.show()
        
    return ax
