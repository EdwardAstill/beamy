from __future__ import annotations
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .member import Member
    from ..beam1d.analysis import LoadedBeam

@dataclass
class MemberResults:
    """Detailed analysis results for an individual member."""
    member: Member
    loaded_beam: LoadedBeam
    
    @property
    def axial(self): return self.loaded_beam.axial()
    @property
    def shear_y(self): return self.loaded_beam.shear("y")
    @property
    def shear_z(self): return self.loaded_beam.shear("z")
    @property
    def bending_y(self): return self.loaded_beam.bending("y")
    @property
    def bending_z(self): return self.loaded_beam.bending("z")
    @property
    def torsion(self): return self.loaded_beam.torsion()
    @property
    def von_mises(self): return self.loaded_beam.von_mises()

