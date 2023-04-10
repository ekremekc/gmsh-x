import gmsh
import os
import sys
import numpy as np
from gmsh_x.mesh_utils import  mirror_mesh_x_axis, fltk_options

mesh_dir = ""
mesh_name = "/Cylinder"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

gmsh.model.add("Geom")
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

path = os.path.dirname(os.path.abspath(__file__))

def copyAndRotate(hole, startAngle, numberOfHoles):
    jumpAngle = 2*np.pi/numberOfHoles
    angle = start_angle
    for i in range(1, numberOfHoles):
        copiedHole = gmsh.model.occ.copy(hole)
        angle += jumpAngle
        print(angle)
        gmsh.model.occ.rotate(copiedHole, 0, 0, 0, 1, 0, 0, angle)
        # gmsh.model.occ.synchronize()

start_angle = 0
N = 40

cylinder = gmsh.model.occ.importShapes('Cylinder.step')
gmsh.model.occ.synchronize()

copyAndRotate(cylinder, start_angle, N)

# copiedHole = gmsh.model.occ.copy(cylinder)
# gmsh.model.occ.rotate(copiedHole, 0, 0, 0, 1, 0, 0, np.pi)
gmsh.model.occ.synchronize()



lc = 0.002
gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10) 
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
gmsh.model.mesh.generate(3)

gmsh.write("{}.msh".format(path+mesh_dir+mesh_name))

if '-nopopup' not in sys.argv:
    fltk_options()
    gmsh.fltk.run()

gmsh.finalize()