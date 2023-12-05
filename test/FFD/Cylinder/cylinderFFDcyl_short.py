from gmsh_x.ffd import FFDCylindrical, deformCylindicalFFD, getMeshdata
from gmsh_x.io_utils import write_xdmf_mesh
import numpy as np
import gmsh
import sys
import os 

if not os.path.exists('MeshDir'):
    os.makedirs('MeshDir')

path = os.path.dirname(os.path.abspath(__file__))
mesh_dir = "/MeshDir"
mesh_name = "/Cylinder"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add(__name__)

R = 0.1

scaler = 0.2
L_total = 1. * scaler

gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L_total, R, tag=1)
gmsh.model.occ.synchronize()

# Physical tags
surfaces = gmsh.model.occ.getEntities(dim=2)

for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

gmsh.model.addPhysicalGroup(3, [1], tag=1) # Geometry tag 

lc = 0.01 #0.005 

gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.Optimize", 1)

gmsh.model.mesh.generate(3)

# Retrieve mesh data before it disappears
mesh_data = getMeshdata(gmsh.model)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("{}.msh".format(path+mesh_dir+mesh_name))

write_xdmf_mesh(path+mesh_dir+mesh_name,dimension=3)

### Introducing FFD 
l, m, n = 2, 4, 5
CylindricalLattice = FFDCylindrical(gmsh.model, l, m , n, 3, tag=-1, includeBoundary=True, parametric=False) 

# Radial Deformation
for i in range(m):
    CylindricalLattice.Pr[1, i, 0] -= 0.02
    CylindricalLattice.Pr[1, i, 1] += 0.02
    CylindricalLattice.Pr[1, i, 2] -= 0.2
    CylindricalLattice.Pr[1, i, 3] += 0.02
    CylindricalLattice.Pr[1, i, 4] -= 0.02

# Angular Deformation
# for i in range(m):
#     CylindricalLattice.Pphi[1, i, 2] += np.pi

CylindricalLattice.write_ffd_points("MeshDir/FFD")

# Start deformation
gmsh.model.add('deformedModel')

gmsh.model = deformCylindicalFFD(gmsh.model, mesh_data, CylindricalLattice)

# Physical tags - using original tags
for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

gmsh.model.addPhysicalGroup(3, [1], tag=1) # Geometry tag 

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/CylinderDeformed"

gmsh.write("{}.msh".format(mesh_name))

write_xdmf_mesh(mesh_name,dimension=3)