import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def plot_extruded_section(section_points, length):
    """
    Plots a 3D extrusion of a 2D section.
    The section is defined in the YZ plane and extruded along the X axis.
    
    Args:
        section_points: List of (y, z) tuples defining the section vertices.
                       Must be a closed loop (first and last point same).
        length: Length of extrusion along the X axis.
    """
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Convert to numpy array for easier handling
    pts = np.array(section_points)
    
    # Validate inputs
    if len(pts) < 3:
        print("Need at least 3 points to define a shape")
        return

    # Vertices at x=0 (Start of beam)
    # Mapping to plot coordinates (User request):
    # Plot X <= Beam Z
    # Plot Y <= Beam X
    # Plot Z <= Beam Y
    x0 = np.zeros(pts.shape[0])
    # (BeamZ, BeamX, BeamY)
    base_verts = np.column_stack((pts[:, 1], x0, pts[:, 0]))

    # Vertices at x=length (End of beam)
    x_end = np.full(pts.shape[0], length)
    end_verts = np.column_stack((pts[:, 1], x_end, pts[:, 0]))

    faces = []

    # 1. Start face (Beam x=0 -> Plot Y=0)
    faces.append(base_verts)

    # 2. End face (Beam x=length -> Plot Y=length)
    faces.append(end_verts)

    # 3. Side faces
    # Connect corresponding points between base and end
    # Iterate through segments. Since points are closed, we have len-1 segments.
    for i in range(len(pts) - 1):
        # Define the quad for the side face
        # Order: Base[i], Base[i+1], End[i+1], End[i]
        side = [
            base_verts[i],
            base_verts[i+1],
            end_verts[i+1],
            end_verts[i]
        ]
        faces.append(side)

    # Create the collection
    mesh = Poly3DCollection(faces, alpha=0.5, edgecolor='k')
    
    # Optional: Set face colors
    mesh.set_facecolor('cyan')

    ax.add_collection3d(mesh)

    # Set plot limits
    # We need to manually set limits because add_collection3d doesn't auto-scale
    all_coords = np.vstack((base_verts, end_verts))
    
    x_min, x_max = all_coords[:, 0].min(), all_coords[:, 0].max()
    y_min, y_max = all_coords[:, 1].min(), all_coords[:, 1].max()
    z_min, z_max = all_coords[:, 2].min(), all_coords[:, 2].max()
    
    # Add some padding
    max_range = np.array([x_max-x_min, y_max-y_min, z_max-z_min]).max() / 2.0
    mid_x = (x_max+x_min) * 0.5
    mid_y = (y_max+y_min) * 0.5
    mid_z = (z_max+z_min) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('Z (Transverse Horizontal)')
    ax.set_ylabel('X (Longitudinal)')
    ax.set_zlabel('Y (Vertical)')
    
    # Set initial view to match the user's request better if needed
    # But the coordinate swap should do the heavy lifting.
    ax.view_init(elev=20, azim=-60) # Default is often 30, -60. 
    
    plt.title(f"Extruded Section (L={length})")
    plt.show()

if __name__ == "__main__":
    # User defined triangle points (y, z)
    points = ((0,0), (2,0), (1,1), (0,0))
    
    # Extrude along X axis
    plot_extruded_section(points, length=5)
