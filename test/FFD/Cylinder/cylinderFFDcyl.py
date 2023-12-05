from math import pi, cos, sin
import gmsh
import os 
import sys
import numpy as np
from math import comb
from gmsh_x.io_utils import write_xdmf_mesh
from gmsh_x.mesh_utils import cyl2cart, cart2cyl

if not os.path.exists('MeshDir'):
    os.makedirs('MeshDir')

path = os.path.dirname(os.path.abspath(__file__))
mesh_dir = "/MeshDir"
mesh_name = "/Cylinder"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add(__name__)

R = 0.1# 0.047/2

scaler = 0.05
L_total = 1. * scaler

gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L_total, R, tag=1)
gmsh.model.occ.synchronize()

# Physical tags
surfaces = gmsh.model.occ.getEntities(dim=2)

for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

gmsh.model.addPhysicalGroup(3, [1], tag=1) # Geometry tag 

lc = 0.1 #0.005 

gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.option.setNumber("Mesh.SaveAll", 0)

gmsh.model.mesh.generate(3)

mesh_data = {}
for e in gmsh.model.getEntities():
    mesh_data[e] = (gmsh.model.getBoundary([e]),
            gmsh.model.mesh.getNodes(e[0], e[1]),
            gmsh.model.mesh.getElements(e[0], e[1]))

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.write("{}.msh".format(path+mesh_dir+mesh_name))

# gmsh.finalize()

write_xdmf_mesh(path+mesh_dir+mesh_name,dimension=3)

### Introducing FFD 

l = 2
m = 4
n = 5
Pr = np.zeros((l,m,n))
Pphi = np.zeros((l,m,n))
Pz = np.zeros((l,m,n))

nodes, coords, param = gmsh.model.mesh.getNodes(3, -1, True, True)
xs = coords[0::3]
ys = coords[1::3]
zs = coords[2::3]

rhos, phis, zetas = cart2cyl(xs, ys, zs)
 
dr = max(rhos)-min(rhos)
dphi = 2*np.pi
dz = max(zetas)-min(zetas)

for i in range(l):
    for j in range(m):
        for k in range(n):
            Pr[i, j, k] = min(rhos)  + dr * i / (l - 1)
            Pphi[i, j, k] = min(phis)  + dphi * j / (m -1)
            Pz[i, j, k] = min(zetas)  + dz * k / (n -1)

P0 = np.array([Pr[0, 0, 0], Pphi[0, 0, 0], Pz[0, 0, 0]])

# Deformation
# for i in range(m):
#     Pr[1, i, 0] -= 0.02
#     Pr[1, i, 1] += 0.02
#     Pr[1, i, 2] -= 0.2
#     Pr[1, i, 3] += 0.02
#     Pr[1, i, 4] -= 0.02
    # Pphi[1, i, 2] += np.pi

edit_double = 0.1
edit_single = 0.1


# Pr[1, 2] += edit_single
# Pr[1, 0] -= edit_double
# Pr[1, -1] -= edit_double
# Py[0, 2] = Py[0, 2] + 0.1
# Px[2, 1] = Px[2, 1] - 0.3

from pyevtk.hl import pointsToVTK

r_ffd = Pr.flatten()
phi_ffd = Pphi.flatten()
z_ffd = Pz.flatten()
x_ffd, y_ffd, z_ffd = cyl2cart(r_ffd, phi_ffd, z_ffd)
pointsToVTK("MeshDir/FFD", x_ffd, y_ffd, z_ffd)

# STU

def calcSTU(coords, P0, dr, dphi, dz):
    """
    Calc STU coordinates
    """
    xs = coords[0::3]
    ys = coords[1::3]
    zs = coords[2::3]

    rhos, phis, zetas = cart2cyl(xs, ys, zs)

    s = (rhos - P0[0])/dr
    t = (phis - P0[1])/dphi
    u = (zetas - P0[2])/dz

    return s,t,u

gmsh.model.add('new_model')

for e in mesh_data:
    # print(mesh_data[e][1][1])
    if len(mesh_data[e][1][1])==3:
        old_coord = mesh_data[e][1][1]

        s,t,u = calcSTU(old_coord,P0, dr, dphi, dz)
        # print(s,t,u)
        Xdef = np.zeros(3)
        for i in range(l):
            for j in range(m):
                for k in range(n):
                    Xdef += comb(l-1,i)*np.power(1-s, l-1-i)*np.power(s,i) * \
                            comb(m-1,j)*np.power(1-t, m-1-j)*np.power(t,j) * \
                            comb(n-1,k)*np.power(1-u, n-1-k)*np.power(u,k) * \
                            np.asarray([Pr[i,j,k], Pphi[i,j,k], Pz[i,j,k]])
        # print(Xdef)
        # Xdef_3d = np.array([Xdef[0],Xdef[1],0])
        Xdef_3d_cart = Xdef.copy()
        Xdef_3d_cart[0], Xdef_3d_cart[1],Xdef_3d_cart[2] = cyl2cart(Xdef[0], Xdef[1], Xdef[2]) 
        new_coord = Xdef_3d_cart.flatten()

    else:
        old_coords = mesh_data[e][1][1]

        s,t,u = calcSTU(old_coords,P0, dr, dphi, dz)
        Xdef = np.zeros((int(len(old_coords)/3),3))
        for point, param_s in enumerate(s):
            for i in range(l):
                for j in range(m):
                    for k in range(n):
                        Xdef[point] +=  comb(l-1,i)*np.power(1-s[point], l-1-i)*np.power(s[point],i) * \
                                        comb(m-1,j)*np.power(1-t[point], m-1-j)*np.power(t[point],j) * \
                                        comb(n-1,k)*np.power(1-u[point], n-1-k)*np.power(u[point],k) * \
                                        np.asarray([Pr[i,j,k], Pphi[i,j,k], Pz[i,j,k]])

        
        Xdef_3d_cart = Xdef.copy()
        Xdef_3d_cart[:,0], Xdef_3d_cart[:,1],Xdef_3d_cart[:,2] = cyl2cart(Xdef[:,0], Xdef[:,1], Xdef[:,2]) 

        new_coord = Xdef_3d_cart.flatten()
        print(Xdef_3d_cart)
        
    gmsh.model.addDiscreteEntity(e[0], e[1], [b[1] for b in mesh_data[e][0]])
    gmsh.model.mesh.addNodes(e[0], e[1], mesh_data[e][1][0], new_coord)
    gmsh.model.mesh.addElements(e[0], e[1], mesh_data[e][2][0], mesh_data[e][2][1], mesh_data[e][2][2])

# Physical tags
for surface in surfaces:
    gmsh.model.addPhysicalGroup(2, [surface[1]])

gmsh.model.addPhysicalGroup(3, [1], tag=1) # Geometry tag 

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

mesh_name = "MeshDir/deformedCylinder"

gmsh.write("{}.msh".format(mesh_name))

write_xdmf_mesh(mesh_name,dimension=3)