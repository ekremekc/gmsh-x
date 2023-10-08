import gmsh
import math
import os
import sys
import numpy as np

gmsh.initialize()

gmsh.model.add("cylinder")

z_perforated_plate = 0
l_perforated_plate = 0.1
r_perforated_plate = 0.1
r_center = 0 

cylinder = gmsh.model.occ.addCylinder(r_center, 0, z_perforated_plate ,0 , 0, l_perforated_plate, r_perforated_plate)
gmsh.model.occ.synchronize()

xmin = 0
ymin = 0
zmin = 0

L = 0.3
H = 0.3

# Create the first cutting plane:
slicer1 = []
slicer1.append((2, gmsh.model.occ.addRectangle(xmin, ymin, zmin, L, H)))
gmsh.model.occ.rotate([slicer1[0]], xmin, ymin, zmin, 0, 1, 0, -math.pi/2)

slicer2 = []
slicer2.append((2, gmsh.model.occ.addRectangle(xmin, ymin, zmin, L, H)))
gmsh.model.occ.rotate([slicer2[0]], xmin, ymin, zmin, 0, 1, 0, -math.pi/2)
gmsh.model.occ.rotate([slicer2[0]], xmin, ymin, zmin, 0, 0, 1, -math.pi/8)

# print(cylinder)
gmsh.model.occ.fragment([(3,cylinder)], slicer1)
gmsh.model.occ.fragment([(3,cylinder)], slicer2)
gmsh.model.occ.remove(gmsh.model.occ.getEntities(2), True)

gmsh.model.occ.synchronize()

print(gmsh.model.occ.getEntities(3))

gmsh.model.occ.remove([(3,1)], True)
gmsh.model.occ.synchronize()

gmsh.option.setNumber('Mesh.MeshSizeMax', 0.005)

gmsh.model.mesh.generate(3)

from gmsh_x.mesh_utils import getMeshData, rotate_z

m = getMeshData(gmsh.model)

angles = np.linspace(360/16, 15*360/16, 15)

for angle in angles:
    print(angle)
    max_node_tag = gmsh.model.mesh.getMaxNodeTag() 
    max_element_tag = gmsh.model.mesh.getMaxElementTag() 
    rotate_z(m, max_node_tag, max_node_tag, max_element_tag, angle)

gmsh.model.occ.synchronize()

gmsh.model.mesh.removeDuplicateNodes()

gmsh.model.occ.synchronize()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()