from typing import Optional
import matplotlib.pyplot as plt
import numpy as np

from .analysis import LoadedBeam

def plot_analysis_results(loaded_beam: LoadedBeam, save_path: Optional[str] = None, show: bool = True, points: int = 100, units: Optional[dict[str, str]] = None):
    """
    Plots the analysis results (Shear, Moment, Deflection, Axial/Torsion) on 2D line graphs.
    
    Args:
        loaded_beam: The analyzed LoadedBeam object.
        save_path: Optional path to save the figure.
        show: Whether to display the plot.
        points: Number of points to sample along the beam.
        units: Optional dictionary specifying units for labels (keys: 'length', 'force', 'moment', 'deflection').
    """
    
    # Setup the figure with 2x2 subplots
    # Changed figsize to be wider (16, 8) instead of (14, 10)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 8))
    # fig.suptitle(f"Beam Analysis Results (L = {loaded_beam.beam.L})", fontsize=16)

    # Process units
    units = units or {}
    u_len = f" ({units['length']})" if 'length' in units else ""
    u_frc = f" ({units['force']})" if 'force' in units else ""
    u_mom = f" ({units['moment']})" if 'moment' in units else ""
    u_def = f" ({units['deflection']})" if 'deflection' in units else ""
    
    # -----------------------------------------------------
    # 1. Shear Force Diagram (Fy, Fz)
    # -----------------------------------------------------
    res_shear_y = loaded_beam.shear("y", points=points)
    res_shear_z = loaded_beam.shear("z", points=points)
    
    # Simplified labels: "y-axis" and "z-axis"
    ax1.plot(res_shear_y.action._x, res_shear_y.action._values, label="y-axis", color="blue")
    ax1.plot(res_shear_z.action._x, res_shear_z.action._values, label="z-axis", color="green", linestyle="--")
    
    ax1.set_title("Shear Force")
    ax1.set_xlabel(f"Position{u_len}") 
    ax1.set_ylabel(f"Force{u_frc}")     
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend()
    ax1.axhline(0, color='black', linewidth=0.5)

    # -----------------------------------------------------
    # 2. Bending Moment Diagram (My, Mz)
    # -----------------------------------------------------
    res_bend_y = loaded_beam.bending("y", points=points) # Moments about Z (Mz) caused by Fy
    res_bend_z = loaded_beam.bending("z", points=points) # Moments about Y (My) caused by Fz
    
    # Note: bending("y") -> returns Mz (moment about z axis due to y-loads)
    #       bending("z") -> returns My (moment about y axis due to z-loads)
    # But user asked for simplified legend. We should probably label them by the axis they are ABOUT, or the axis of bending?
    # Usually "y-axis" bending means bending about Y or bending in Y-plane?
    # User said: "so shear Fy (vertical just becomes y-axis)"
    # Shear Fy is force in Y direction.
    # Moment Mz is moment about Z direction (associated with Fy).
    # If I label Mz as "z-axis", it matches the moment vector direction.
    
    ax2.plot(res_bend_y.action._x, res_bend_y.action._values, label="z-axis", color="blue")
    ax2.plot(res_bend_z.action._x, res_bend_z.action._values, label="y-axis", color="green", linestyle="--")
    
    ax2.set_title("Bending Moment")
    ax2.set_xlabel(f"Position{u_len}")
    ax2.set_ylabel(f"Moment{u_mom}")
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()
    ax2.axhline(0, color='black', linewidth=0.5)

    # -----------------------------------------------------
    # 3. Deflection (dy, dz)
    # -----------------------------------------------------
    disp_y = loaded_beam.deflection("y", points=points) # Deflection in y direction
    disp_z = loaded_beam.deflection("z", points=points) # Deflection in z direction
    
    ax3.plot(disp_y._x, disp_y._values, label="y-axis", color="blue")
    ax3.plot(disp_z._x, disp_z._values, label="z-axis", color="green", linestyle="--")
    
    ax3.set_title("Deflection")
    ax3.set_xlabel(f"Position{u_len}")
    ax3.set_ylabel(f"Displacement{u_def}")
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend()
    ax3.axhline(0, color='black', linewidth=0.5)
    
    # -----------------------------------------------------
    # 4. Axial Force and Torsion
    # -----------------------------------------------------
    res_axial = loaded_beam.axial(points=points)
    res_torsion = loaded_beam.torsion(points=points)
    
    # Plot Axial on left axis
    line1 = ax4.plot(res_axial.action._x, res_axial.action._values, label="Axial Force", color="red")
    ax4.set_ylabel(f"Axial Force{u_frc}", color="red") 
    ax4.tick_params(axis='y', labelcolor="red")
    
    # Plot Torsion on right axis
    ax4_right = ax4.twinx()
    line2 = ax4_right.plot(res_torsion.action._x, res_torsion.action._values, label="Torsion", color="purple", linestyle="--")
    ax4_right.set_ylabel(f"Torsion Moment{u_mom}", color="purple")
    ax4_right.tick_params(axis='y', labelcolor="purple")
    
    ax4.set_title("Axial Force & Torsion")
    ax4.set_xlabel(f"Position{u_len}") 
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.axhline(0, color='black', linewidth=0.5)
    
    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='upper right')

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    
    if show:
        plt.show()
