from typing import Optional, Literal
import matplotlib.pyplot as plt
import numpy as np

from sectiony.stress import Stress, StressType
from .analysis import LoadedBeam

class StressPlotter:
    """
    Plots stress distributions on the beam's cross-section at a specific location along the beam.
    """
    def __init__(self, loaded_beam: LoadedBeam):
        self.loaded_beam = loaded_beam
        self.beam = loaded_beam.beam

    def _get_internal_forces_at(self, x_pos: float) -> dict[str, float]:
        """
        Interpolate internal forces at position x_pos.
        
        Returns:
            Dictionary with keys: N, Vy, Vz, Mx, My, Mz
        """
        # Get analysis results
        # These methods return (x_coords, values)
        # We need to interpolate at x_pos
        
        # Axial Force (N) = Fx
        res_axial = self.loaded_beam.axial(points=100) # Use default points or analyze on fly?
        # Ideally we want exact values. LoadedBeam methods analyze at n points.
        # However, LoadedBeam.axial() returns AnalysisResult which contains .action (forces)
        # We can use the interpolation functions used inside AnalysisResult, but they are private or internal.
        # Actually LoadedBeam methods run the analysis if not cached, but they return arrays.
        # A better way is to use the underlying interpolation functions from analysis.py if possible,
        # but LoadedBeam abstracts that.
        # Let's use the .at(x) method if AnalysisResult had one, but it doesn't seem to.
        # It returns a Result object which has ._x and ._values. We can interpolate from that.
        
        # Helper to interpolate from Result object
        def get_val(analysis_res, x):
            return np.interp(x, analysis_res.action._x, analysis_res.action._values)

        N = get_val(self.loaded_beam.axial(), x_pos)
        
        # Shear Forces (Vy, Vz)
        # Note: Beamy's sign conventions might need mapping to Sectiony's
        # Beamy: Y is vertical, Z is horizontal? 
        # Check coordinates:
        # Beamy 1D: x is along beam. 
        # Beamy Section: y, z are cross section coordinates.
        # Sectiony: y, z are cross section coordinates.
        # Assuming they align.
        
        Vy = get_val(self.loaded_beam.shear("y"), x_pos)
        Vz = get_val(self.loaded_beam.shear("z"), x_pos)
        
        # Torsion (Mx)
        Mx = get_val(self.loaded_beam.torsion(), x_pos)
        
        # Bending Moments (My, Mz)
        # Beamy bending("y") corresponds to transverse loads in Y, which cause bending moments about Z (Mz).
        # Beamy bending("z") corresponds to transverse loads in Z, which cause bending moments about Y (My).
        
        # Map Beamy moments to Sectiony moments:
        # Sectiony My needs Moment about Y -> Beamy bending("z")
        # Sectiony Mz needs Moment about Z -> Beamy bending("y")
        
        My = get_val(self.loaded_beam.bending("z"), x_pos)
        Mz = get_val(self.loaded_beam.bending("y"), x_pos)
        
        return {
            "N": N,
            "Vy": Vy,
            "Vz": Vz,
            "Mx": Mx,
            "My": My,
            "Mz": Mz
        }

    def plot_stress_at(
        self, 
        x_pos: float, 
        stress_type: StressType = "von_mises", 
        ax: Optional[plt.Axes] = None, 
        show: bool = True,
        cmap: str = "viridis",
        title: Optional[str] = None
    ) -> Optional[plt.Axes]:
        """
        Plot stress distribution on the section at a specific beam location x.
        
        Args:
            x_pos: Position along the beam length [0, L].
            stress_type: Type of stress to plot (e.g., 'von_mises', 'sigma_bending', etc.).
            ax: Matplotlib axes to plot on.
            show: Whether to show the plot immediately.
            cmap: Colormap to use.
            title: Optional title override.
        """
        if not (0 <= x_pos <= self.beam.L):
            print(f"Warning: Position {x_pos} is outside beam length [0, {self.beam.L}]")
            # We can still try to plot, or clamp, or just proceed (interp will extrapolate or clamp)
        
        forces = self._get_internal_forces_at(x_pos)
        
        # Create Sectiony Stress object
        stress = Stress(
            section=self.beam.section,
            N=forces["N"],
            Vy=forces["Vy"],
            Vz=forces["Vz"],
            Mx=forces["Mx"],
            My=forces["My"],
            Mz=forces["Mz"]
        )
        
        # Plot
        ax = stress.plot(stress_type=stress_type, ax=ax, show=False, cmap=cmap)
        
        if ax:
            # Add context to title
            if title is None:
                current_title = ax.get_title()
                ax.set_title(f"{current_title}\n@ x={x_pos:.3f}")
            else:
                ax.set_title(title)
                
            if show:
                plt.show()
                
        return ax

