from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional

def _polygon_area_signed(points: List[Tuple[float, float]]) -> float:
    area = 0.0
    for i in range(len(points)):
        j = (i + 1) % len(points)
        area += points[i][0] * points[j][1]
        area -= points[j][0] * points[i][1]
    return 0.5 * area

def _polygon_centroid(points: List[Tuple[float, float]]) -> Tuple[float, float]:
    area = _polygon_area_signed(points)
    if abs(area) < 1e-9: return (0.0, 0.0)
    cx = 0.0
    cy = 0.0
    for i in range(len(points)):
        j = (i + 1) % len(points)
        factor = points[i][0] * points[j][1] - points[j][0] * points[i][1]
        cx += (points[i][0] + points[j][0]) * factor
        cy += (points[i][1] + points[j][1]) * factor
    return (cx / (6 * area), cy / (6 * area))

def _polygon_properties(points: List[Tuple[float, float]]) -> Tuple[float, float, float, float, float]:
    A_signed = _polygon_area_signed(points)
    A = abs(A_signed)
    if A < 1e-9:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    
    Cy, Cz = _polygon_centroid(points)
    
    Iyy_origin = 0.0 # integral z^2 dA
    Izz_origin = 0.0 # integral y^2 dA
    
    for i in range(len(points)):
        j = (i + 1) % len(points)
        y1, z1 = points[i]
        y2, z2 = points[j]
        cross = y1 * z2 - y2 * z1
        
        Izz_origin += (y1**2 + y1*y2 + y2**2) * cross
        Iyy_origin += (z1**2 + z1*z2 + z2**2) * cross
        
    Izz_origin = abs(Izz_origin) / 12.0
    Iyy_origin = abs(Iyy_origin) / 12.0
    
    # Parallel axis theorem to centroid
    Iyy = Iyy_origin - A * Cz**2
    Izz = Izz_origin - A * Cy**2
    
    ys = [p[0] for p in points]
    zs = [p[1] for p in points]
    
    y_max_dist = max(ys) - Cy
    z_max_dist = max(zs) - Cz
    
    return A, Iyy, Izz, y_max_dist, z_max_dist

def _clip_polygon(subject: List[Tuple[float, float]], clipper: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    # Sutherland-Hodgman clipping
    def inside(p, cp1, cp2):
        return (cp2[0]-cp1[0])*(p[1]-cp1[1]) - (cp2[1]-cp1[1])*(p[0]-cp1[0]) >= 0
    
    def compute_intersection(cp1, cp2, s, e):
        dc = [cp1[0] - cp2[0], cp1[1] - cp2[1]]
        dp = [s[0] - e[0], s[1] - e[1]]
        n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
        n2 = s[0] * e[1] - s[1] * e[0]
        det = dc[0] * dp[1] - dc[1] * dp[0]
        if abs(det) < 1e-9: return e
        n3 = 1.0 / det
        return ((n1 * dp[0] - n2 * dc[0]) * n3, (n1 * dp[1] - n2 * dc[1]) * n3)
    
    output = list(subject)
    if _polygon_area_signed(clipper) < 0:
        clipper_oriented = list(reversed(clipper))
    else:
        clipper_oriented = list(clipper)
        
    cp1 = clipper_oriented[-1]
    for cp2 in clipper_oriented:
        input_list = output
        output = []
        if not input_list:
            break
        s = input_list[-1]
        for e in input_list:
            if inside(e, cp1, cp2):
                if not inside(s, cp1, cp2):
                    output.append(compute_intersection(cp1, cp2, s, e))
                output.append(e)
            elif inside(s, cp1, cp2):
                output.append(compute_intersection(cp1, cp2, s, e))
            s = e
        cp1 = cp2
    return output

@dataclass
class Shape:
    points: List[Tuple[float, float]]
    hollow: bool = False

@dataclass
class Geometry:
    shapes: List[Shape]

    def reduce_hollows(self) -> List[Shape]:
        solids = [s for s in self.shapes if not s.hollow]
        hollows = [s for s in self.shapes if s.hollow]
        
        if not hollows:
            return solids
            
        reduced_shapes = list(solids)
        for h in hollows:
            for s in solids:
                clipped_points = _clip_polygon(h.points, s.points)
                if len(clipped_points) >= 3 and abs(_polygon_area_signed(clipped_points)) > 1e-9:
                    reduced_shapes.append(Shape(points=clipped_points, hollow=True))
        return reduced_shapes


@dataclass
class Section:
    """
    Cross-section properties in local coordinates.
    Can be initialized with properties directly or with a geometry object.
    """
    name: str
    A: Optional[float] = None
    Iy: Optional[float] = None
    Iz: Optional[float] = None
    J: Optional[float] = None
    y_max: Optional[float] = None
    z_max: Optional[float] = None
    geometry: Optional[Geometry] = None

    def plot(self, ax=None, show=True):
        """Plot the section geometry."""
        from .plotter import plot_section
        return plot_section(self, ax=ax, show=show)

    def __post_init__(self):
        if self.geometry and self.A is None:
            self._calculate_properties_from_geometry()
        
        if self.A is None:
            raise ValueError("Section properties must be provided if geometry is not.")

    def _calculate_properties_from_geometry(self):
        shapes = self.geometry.reduce_hollows()
        total_A = 0.0
        total_Ay = 0.0
        total_Az = 0.0
        
        for s in shapes:
            A, _, _, _, _ = _polygon_properties(s.points)
            sign = -1.0 if s.hollow else 1.0
            cx, cz = _polygon_centroid(s.points)
            total_A += sign * A
            total_Ay += sign * A * cx
            total_Az += sign * A * cz
            
        if total_A == 0:
            self.A = 0.0
            self.Iy = 0.0
            self.Iz = 0.0
            self.J = 0.0
            self.y_max = 0.0
            self.z_max = 0.0
            return
            
        Cy_global = total_Ay / total_A
        Cz_global = total_Az / total_A
        
        total_Iyy = 0.0
        total_Izz = 0.0
        max_y_dist = 0.0
        max_z_dist = 0.0
        
        for s in shapes:
            A, Iyy_local, Izz_local, _, _ = _polygon_properties(s.points)
            sign = -1.0 if s.hollow else 1.0
            cx, cz = _polygon_centroid(s.points)
            
            dy = cx - Cy_global
            dz = cz - Cz_global
            
            total_Iyy += sign * (Iyy_local + A * dz**2)
            total_Izz += sign * (Izz_local + A * dy**2)
            
            if not s.hollow:
                ys = [p[0] for p in s.points]
                zs = [p[1] for p in s.points]
                max_y_dist = max(max_y_dist, max(ys) - Cy_global)
                max_z_dist = max(max_z_dist, max(zs) - Cz_global)
                
        self.A = total_A
        self.Iy = total_Iyy
        self.Iz = total_Izz
        self.J = total_Iyy + total_Izz # Approx polar moment
        self.y_max = max_y_dist
        self.z_max = max_z_dist

