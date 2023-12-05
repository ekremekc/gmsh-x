import gmsh
import sys
import numpy as np
from math import comb
from gmsh_x.io_utils import write_xdmf_mesh
from gmsh_x.mesh_utils import cyl2cart, cart2cyl

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

mesh_name = "MeshDir/Circle"

gmsh.write("{}.msh".format(mesh_name))

write_xdmf_mesh(mesh_name,dimension=2)

### Introducing FFD 

l = 2
m = 4
Pr = np.zeros((l,m))
Pphi = np.zeros((l,m))

nodes, coords, param = gmsh.model.mesh.getNodes(2, s1, True, True)
xs = coords[0::3]
ys = coords[1::3]
zs = coords[1::3]

rhos, phis, zetas = cart2cyl(xs, ys, zs)
 
dr = max(rhos)-min(rhos)
dphi = 2*np.pi
for i in range(l):
    for j in range(m):
        Pr[i, j] = min(rhos)  + dr * i / (l - 1)
        Pphi[i, j] = min(phis)  + dphi * j / (m -1)

P0 = np.array([Pr[0, 0], Pphi[0, 0]])

# Deformation
# for i in range(m):
#     Pr[1, i] += 0.00
#     Pphi[0, i] += np.pi/2

edit_double = 0.1
edit_single = 0.1

Pr[1, 0] += edit_double
Pr[1, -1] += edit_double

from pyevtk.hl import pointsToVTK

r_ffd = Pr.flatten()
phi_ffd = Pphi.flatten()
z_ffd = np.zeros_like(phi_ffd)
x_ffd, y_ffd, z_ffd = cyl2cart(r_ffd, phi_ffd, z_ffd)
pointsToVTK("MeshDir/FFD", x_ffd, y_ffd, z_ffd)

# STU

def calcSTU(coords, P0, dr, dphi):
    """
    Calc STU coordinates
    """
    xs = coords[0::3]
    ys = coords[1::3]
    zs = coords[2::3]

    rhos, phis, zetas = cart2cyl(xs, ys, zs)

    s = (rhos - P0[0])/dr
    t = (phis - P0[1])/dphi

    return s,t

gmsh.model.add('new_model')

for e in mesh_data:
    # print(mesh_data[e][1][1])
    if len(mesh_data[e][1][1])==3:
        old_coord = mesh_data[e][1][1]
        s,t = calcSTU(old_coord,P0, dr, dphi)
        Xdef = np.zeros(2)
        for i in range(l):
            for j in range(m):
                Xdef +=  comb(l-1,i)*np.power(1-s, l-1-i)*np.power(s,i) * \
                                comb(m-1,j)*np.power(1-t, m-1-j)*np.power(t,j) * \
                                    np.asarray([Pr[i,j], Pphi[i,j]])

        Xdef_3d = np.array([Xdef[0],Xdef[1],0])
        Xdef_3d_cart = Xdef_3d.copy()
        Xdef_3d_cart[0], Xdef_3d_cart[1],Xdef_3d_cart[2] = cyl2cart(Xdef_3d[0], Xdef_3d[1], Xdef_3d[2]) 
        new_coord = Xdef_3d_cart.flatten()

    else:
        old_coords = mesh_data[e][1][1]

        s,t = calcSTU(old_coords,P0, dr, dphi)
        Xdef = np.zeros((int(len(old_coords)/3),2))
        for point, param_s in enumerate(s):
            for i in range(l):
                for j in range(m):
                    Xdef[point] +=  comb(l-1,i)*np.power(1-s[point], l-1-i)*np.power(s[point],i) * \
                                    comb(m-1,j)*np.power(1-t[point], m-1-j)*np.power(t[point],j) * \
                                        np.asarray([Pr[i,j], Pphi[i,j]])

        
        z_axis = np.zeros((len(Xdef),1))
        Xdef_3d = np.hstack((Xdef,z_axis))

        Xdef_3d_cart = Xdef_3d.copy()
        Xdef_3d_cart[:,0], Xdef_3d_cart[:,1],Xdef_3d_cart[:,2] = cyl2cart(Xdef_3d[:,0], Xdef_3d[:,1], Xdef_3d[:,2]) 

        new_coord = Xdef_3d_cart.flatten()
        # print(Xdef_3d_cart)
        
    gmsh.model.addDiscreteEntity(e[0], e[1], [b[1] for b in mesh_data[e][0]])
    gmsh.model.mesh.addNodes(e[0], e[1], mesh_data[e][1][0], new_coord)
    gmsh.model.mesh.addElements(e[0], e[1], mesh_data[e][2][0], mesh_data[e][2][1], mesh_data[e][2][2])

gmsh.model.addPhysicalGroup(1,[1],1)
gmsh.model.addPhysicalGroup(2,[1],1)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/deformedCircle"

gmsh.write("{}.msh".format(mesh_name))

write_xdmf_mesh(mesh_name,dimension=2)