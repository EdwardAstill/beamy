from typing import Optional
import matplotlib.pyplot as plt
from sectiony import Section
from sectiony.plotter import plot_section as sectiony_plot_section

def plot_section(section: Section, ax: Optional[plt.Axes] = None, show: bool = True) -> Optional[plt.Axes]:
    """
    Plot the cross-section geometry.
    Wrapper around sectiony.plotter.plot_section.
    
    Args:
        section: The section to plot
        ax: Optional matplotlib axes to plot on. If None, creates a new figure.
        show: Whether to call plt.show() at the end.
    
    Returns:
        The axes object used for plotting.
    """
    return sectiony_plot_section(section, ax=ax, show=show)

