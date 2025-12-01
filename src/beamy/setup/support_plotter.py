from typing import TYPE_CHECKING, Optional, List
import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from .beam import Support

def plot_supports(supports: List["Support"], beam_length: float, unit: str = "m", save_path: Optional[str] = None, show: bool = True):
    """
    Plots the beam as a straight line with supports marked as dots and labeled.
    
    Args:
        supports: List of Support objects.
        beam_length: Length of the beam.
        unit: Optional unit string to append to position labels (e.g., "m", "mm").
        save_path: Optional path to save the figure.
        show: Whether to display the plot.
    """
    fig, ax = plt.subplots(figsize=(10, 2))
    
    # Draw beam axis line explicitly (since we hide the spine)
    ax.plot([0, beam_length], [0, 0], color='black', linewidth=2, zorder=1)
    
    # Plot supports
    if supports:
        xs = [s.x for s in supports]
        ys = [0] * len(supports)
        types = [s.type for s in supports]
        
        # Plot markers
        ax.scatter(xs, ys, marker='o', color='red', s=100, zorder=2, label='Supports')
        
        # Add support type labels (above)
        for x, y, t in zip(xs, ys, types):
            ax.annotate(t, (x, y), xytext=(0, 15), textcoords='offset points', 
                        ha='center', va='bottom', rotation=45)
            
            # Add position labels (below)
            label = f"{x} {unit}" if unit else f"{x}"
            ax.annotate(label, (x, y), xytext=(0, -20), textcoords='offset points',
                        ha='center', va='top')

    # Styling
    # Restore wider limits to see labels clearly, but centered on beam
    margin = beam_length * 0.1
    ax.set_xlim(-margin, beam_length + margin)
    ax.set_ylim(-1, 1)
    
    # Hide everything else
    ax.axis('off')
    ax.set_title('Beam Supports')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    
    if show:
        plt.show()
