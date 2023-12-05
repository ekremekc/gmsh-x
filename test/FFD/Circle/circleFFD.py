import gmsh
import sys
import numpy as np
from math import comb

gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("starting_model")

r = 0.2
circle = gmsh.model.occ.addCircle(0, 0, 0, r=r)
cl = gmsh.model.occ.addCurveLoop([circle])
s1 = gmsh.model.occ.addPlaneSurface([cl])

gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeMin",r/20)
gmsh.option.setNumber("Mesh.MeshSizeMax",r/10)
gmsh.model.mesh.generate(2)

mesh_data = {}
for e in gmsh.model.getEntities():
    mesh_data[e] = (gmsh.model.getBoundary([e]),
            gmsh.model.mesh.getNodes(e[0], e[1]),
            gmsh.model.mesh.getElements(e[0], e[1]))

gmsh.model.addPhysicalGroup(1,[1],1)
gmsh.model.addPhysicalGroup(2,[1],1)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/original_circle"

gmsh.write("{}.msh".format(mesh_name))

from helmholtz_x.dolfinx_utils import write_xdmf_mesh
write_xdmf_mesh(mesh_name,dimension=2)

### Introducing FFD 

l = m = 3
Px = np.zeros((l,m))
Py = np.zeros((l,m))

nodes, coords, param = gmsh.model.mesh.getNodes(2, s1, True, True)
xs = coords[0::3]
ys = coords[1::3]

dx = max(xs)-min(xs)
dy = max(ys)-min(ys)

for i in range(l):
    for j in range(m):
        Px[i, j] = min(xs)  + dx * i / (l - 1)
        Py[i, j] = min(ys)  + dy * j / (m - 1)

P0 = np.array([Px[0, 0], Py[0, 0]])

# Deformation
Py[0, 0] = Py[0, 0] - 0.1
Py[0, 2] = Py[0, 2] + 0.1
Px[2, 1] = Px[2, 1] - 0.3

# STU

def calcSTU(coords, P0, dx, dy):
    """
    Calc STU coordinates
    """
    xs = coords[0::3]
    ys = coords[1::3]

    s = (xs - P0[0])/dx
    t = (ys - P0[1])/dy

    return s,t

gmsh.model.add('new_model')

for e in mesh_data:
    # print(mesh_data[e][1][1])
    if len(mesh_data[e][1][1])==3:
        old_coord = mesh_data[e][1][1]
        s = (old_coord[0] - P0[0])/dx
        t = (old_coord[1] - P0[1])/dy
        # print(s,t)
        Xdef = np.zeros(2)
        for i in range(l):
            for j in range(m):
                Xdef +=  comb(l-1,i)*np.power(1-s, l-1-i)*np.power(s,i) * \
                                comb(m-1,j)*np.power(1-t, m-1-j)*np.power(t,j) * \
                                    np.asarray([Px[i,j], Py[i,j]])
        # print(Xdef)
        new_coord = np.array([Xdef[0],Xdef[1],0])
    else:
        old_coords = mesh_data[e][1][1]
        s,t = calcSTU(old_coords,P0, dx, dy)
        print(s)
        old_coords = old_coords.reshape(-1,3)
        Xdef = np.zeros((len(old_coords),2))
        print(old_coords)
        for point, param_s in enumerate(s):
            for i in range(l):
                for j in range(m):
                    Xdef[point] +=  comb(l-1,i)*np.power(1-s[point], l-1-i)*np.power(s[point],i) * \
                                    comb(m-1,j)*np.power(1-t[point], m-1-j)*np.power(t[point],j) * \
                                        np.asarray([Px[i,j], Py[i,j]])
        z_axis = np.zeros((len(Xdef),1))
        Xdef_3d = np.hstack((Xdef,z_axis))
        new_coord = Xdef_3d.flatten()
        print(Xdef_3d)
        
    gmsh.model.addDiscreteEntity(e[0], e[1], [b[1] for b in mesh_data[e][0]])
    gmsh.model.mesh.addNodes(e[0], e[1], mesh_data[e][1][0], new_coord)
    gmsh.model.mesh.addElements(e[0], e[1], mesh_data[e][2][0], mesh_data[e][2][1], mesh_data[e][2][2])

gmsh.model.addPhysicalGroup(1,[1],1)
gmsh.model.addPhysicalGroup(2,[1],1)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/deformed_circle"

gmsh.write("{}.msh".format(mesh_name))

from helmholtz_x.dolfinx_utils import write_xdmf_mesh
write_xdmf_mesh(mesh_name,dimension=2)