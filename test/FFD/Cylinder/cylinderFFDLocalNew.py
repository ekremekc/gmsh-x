from gmsh_x.ffd import FFDCylindrical, getLocalMeshdata, getNonLocalMeshdata, deformCylindicalLocalFFD, getLocalMeshdataNew
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

R = 0.047/2

L_total = 0.5 / 5

gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L_total/2, R)
gmsh.model.occ.synchronize()

gmsh.model.occ.addCylinder(0, 0, L_total/2, 0, 0, L_total/2, R)
gmsh.model.occ.synchronize()

gmsh.model.occ.addCylinder(0, 0, L_total, 0, 0, L_total/2, R)
gmsh.model.occ.synchronize()

gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

# Physical tags
surfaces = gmsh.model.occ.getEntities(dim=2)
for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

volumes = gmsh.model.occ.getEntities(dim=3)
for voltag, volume in enumerate(volumes):
    gmsh.model.addPhysicalGroup(3, [volume[1]], voltag)

lc = 0.01 #0.005 

gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.Optimize", 1)

gmsh.model.mesh.generate(3)

ffd_volume_tag = 1
# Retrieve mesh data before it disappears
# mesh_data_local = getLocalMeshdata(gmsh.model, 3, ffd_volume_tag)
# mesh_data_nonlocal = getNonLocalMeshdata(gmsh.model, 3, ffd_volume_tag)

print("\n\n")

mesh_data_local, mesh_data_nonlocal = getLocalMeshdataNew(3, 2)

print(mesh_data_local)
print(mesh_data_nonlocal)


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("{}.msh".format(path+mesh_dir+mesh_name))

write_xdmf_mesh(path+mesh_dir+mesh_name,dimension=3)

### Introducing FFD 
l, m, n = 2, 4, 5
CylindricalLattice = FFDCylindrical(gmsh.model, l, m , n, 3, tag=ffd_volume_tag, includeBoundary=True, parametric=False) 

# Radial Deformation
for i in range(m):
    CylindricalLattice.Pr[1, i, 2] += 0.01
    # CylindricalLattice.Pr[1, i, 1] += 0.02
    # CylindricalLattice.Pr[1, i, 2] -= 0.02
    # CylindricalLattice.Pr[1, i, 3] += 0.02
    # CylindricalLattice.Pr[1, i, 4] -= 0.02

# Angular Deformation
# for i in range(m):
#     CylindricalLattice.Pphi[1, i, 2] += np.pi

CylindricalLattice.write_ffd_points("MeshDir/LocalFFD")

# Start deformation
gmsh.model.add('deformedModel')

gmsh.model = deformCylindicalLocalFFD(gmsh.model, mesh_data_local, mesh_data_nonlocal, CylindricalLattice)

# Physical tags - using original tags
for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

for voltag, volume in enumerate(volumes):
    gmsh.model.addPhysicalGroup(3, [volume[1]], voltag)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/CylinderDeformed"

gmsh.write("{}.msh".format(mesh_name))

write_xdmf_mesh(mesh_name,dimension=3)