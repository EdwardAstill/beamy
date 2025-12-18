# Empty structure:
# Load block:
# Sling:
# Total:

# Position of padeyes:
# 1:
# 2:
# 3:
# 4:

# Lifting point position:

# Forces
# 1: (...,...,...) 
# 2: (...,...,...)
# 3: (...,...,...)
# 4: (...,...,...)

# Sum of forces in z direction must equal the wight
import math
#Masses (kg):

empty_structure = 1041
load_block = 2674
sling = 13.88 *46
total = empty_structure + load_block + sling
weight_force = total * 9.81  # in Newtons



# Elevated
print("="*40)
print("Elevated:")
print("="*40)

max_pin_distance = (75 + 10)  # mm
print(f"Max pin distance: {max_pin_distance} mm")

force_on_each_col = (load_block * 9.81)/4  # N
print(f"Force on each column: {force_on_each_col} N")

moment_on_each_beam = force_on_each_col * max_pin_distance  # Nmm
print(f"Moment on each beam: {moment_on_each_beam} Nmm")

applied_height = 1062 - 26 + 75  # mm (from bottom)
print(f"Applied height: {applied_height} mm")


# Lifted
print("="*40)
print("Lifted:")
print("="*40)

max_pin_distance = (75 + 10)  # mm
print(f"Max pin distance: {max_pin_distance} mm")

force_on_each_col = (load_block * 9.81)/4  # N
print(f"Force on each column: {force_on_each_col} N")

moment_on_each_beam = force_on_each_col * max_pin_distance  # Nmm
print(f"Moment on each beam: {moment_on_each_beam} Nmm")

applied_height = 875  # mm (from bottom)
print(f"Applied height: {applied_height} mm")

lift_height = 1200  # mm

pad_to_centre = math.sqrt(500**2 + 400**2)

lift_angle = (math.atan2(lift_height, pad_to_centre)) * (180/math.pi)
print(f"Lift angle: {lift_angle} degrees")

sling_and_structure_weight = (empty_structure + sling) * 9.81

print(f"Base force (sling and structure weight): {sling_and_structure_weight} N")









# # Position of padeyes (x,y,z) in mm)
# padeye_1 = (-500,-425,0)
# padeye_2 = (500,-425,0)
# padeye_3 = (500,425,0)
# padeye_4 = (-500,425,0)

# # Lifting point position (x,y,z) in mm
# lifting_point = (0,0,1200)


# # Force direction vectors
# def force_vector(padeye, lifting_point, weight_force):
#     vector = (lifting_point[0] - padeye[0],
#               lifting_point[1] - padeye[1],
#               lifting_point[2] - padeye[2])
#     magnitude = (vector[0]**2 + vector[1]**2 + vector[2]**2) ** 0.5
#     unit_vector = (vector[0]/magnitude, vector[1]/magnitude, vector[2]/magnitude)
#     # scale vector so that their z component is weight force/4
#     scale = (weight_force / 4) / unit_vector[2]
#     return (unit_vector[0]*scale, unit_vector[1]*scale, unit_vector[2]*scale)
# force_1 = force_vector(padeye_1, lifting_point, weight_force)
# force_2 = force_vector(padeye_2, lifting_point, weight_force)
# force_3 = force_vector(padeye_3, lifting_point, weight_force)
# force_4 = force_vector(padeye_4, lifting_point, weight_force)
# for force in [force_1, force_2, force_3, force_4]:
#     print(f"Force: {force}")

# lift_angle = (math.atan2(lifting_point[2], ((lifting_point[0]-padeye_1[0])**2 + (lifting_point[1]-padeye_1[1])**2)**0.5)) * (180/math.pi)

# print(f"Lift angle: {lift_angle} degrees")