import gmsh
import numpy as np
from math import comb
from pyevtk.hl import pointsToVTK

def cart2cyl(x, y, z):
    # cylindrical to Cartesian
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    zeta = z
    return rho, phi, zeta 

class FFDBox:
    def __init__(self, gmsh_model, l, m , n, dim, tag=-1):
        """This class generates a box-shaped FFD Lattice 

        Args:
            gmsh_model (gmsh.model): gmsh model (after meshing) 
            l (int): number of points in the x direction
            m (int): number of points in the y direction
            n (int): number of points in the z direction
            dim (int): dimension of gmsh entity - 2 or 3
            tag (int, optional): tag of the entity, -1 returns all tags. Defaults to -1.
        """
        self.l = l
        self.m = m
        self.n = n

        self.Px = np.zeros((l,m,n))
        self.Py = np.zeros((l,m,n))
        self.Pz = np.zeros((l,m,n))

        nodes, coords, param = gmsh_model.mesh.getNodes(dim, tag, True, True)

        self.base_coords = coords

        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        self.dx = max(xs)-min(xs)
        self.dy = max(ys)-min(ys)
        self.dz = max(zs)-min(zs)

        for i in range(l):
            for j in range(m):
                for k in range(n):
                    self.Px[i, j, k] = min(xs)  + self.dx * i / (l - 1)
                    self.Py[i, j, k] = min(ys)  + self.dy * j / (m - 1)
                    self.Pz[i, j, k] = min(zs)  + self.dz * k / (n - 1)

        self.P0 = np.array([self.Px[0, 0, 0], self.Py[0, 0, 0], self.Pz[0, 0, 0]])

    def write_ffd_points(self, name="MeshDirFFD/FFDPoints"):
        """writes the FFD points as a vtu file (readable using ParaView).
        """
        
        x_ffd = self.Px.flatten()
        y_ffd = self.Py.flatten()
        z_ffd = self.Pz.flatten()
        pointsToVTK(name, x_ffd, y_ffd, z_ffd)
    
    def calcSTU(self, coords):
        """
        Calc STU (parametric) coordinates
        """
        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        s = (xs - self.P0[0])/self.dx
        t = (ys - self.P0[1])/self.dy
        u = (zs - self.P0[2])/self.dz

        return s,t,u

class FFDAngular:
    def __init__(self, gmsh_model, l, m , n, dim, tag=-1, includeBoundary=False, parametric=True):
        """This class generates a angular-shaped FFD Lattice 

        Args:
            gmsh_model (gmsh.model): gmsh model (after meshing) 
            l (int): number of points in the x direction
            m (int): number of points in the y direction
            n (int): number of points in the z direction
            dim (int): dimension of gmsh entity - 2 or 3
            tag (int, optional): tag of the entity, -1 returns all tags. Defaults to -1.
        """
        self.l = l
        self.m = m
        self.n = n

        self.Px = np.zeros((l,m,n))
        self.Py = np.zeros((l,m,n))
        self.Pz = np.zeros((l,m,n))

        nodes, coords, param = gmsh_model.mesh.getNodes(dim, tag, includeBoundary, parametric)

        self.base_coords = coords

        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        rhos, phis, zetas = cart2cyl(xs, ys, zs)

        self.dr = max(rhos)-min(rhos)
        self.dphi = max(phis)-min(phis)
        self.dz = max(zs)-min(zs)

        for i in range(l):
            for j in range(m):
                for k in range(n):
                    rho = min(rhos)  + self.dr * i / (l - 1)
                    phi = min(phis)  + self.dphi * j / (m - 1)
                    x = rho * np.cos(phi)
                    y = rho * np.sin(phi)
                    self.Px[i, j, k] = x
                    self.Py[i, j, k] = y
                    self.Pz[i, j, k] = min(zs)  + self.dz * k / (n - 1)

        self.P0 = np.array([self.Px[0, 0, 0], self.Py[0, 0, 0], self.Pz[0, 0, 0]])
        
        self.Pr, self.Pphi, self.Pz = cart2cyl(self.Px, self.Py, self.Pz)

    def write_ffd_points(self, name):
        """writes the FFD points as a vtu file (readable using ParaView).
        """
        
        x_ffd = self.Px.flatten()
        y_ffd = self.Py.flatten()
        z_ffd = self.Pz.flatten()
        pointsToVTK(name, x_ffd, y_ffd, z_ffd)
        print("FFD points are saved.")
    
    def calcSTU(self, coords):
        """
        Calc STU (parametric) coordinates
        """
        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        rhos, phis, zetas = cart2cyl(xs, ys, zs)
        rho0, phi0, z0 = cart2cyl(self.P0[0], self.P0[1], self.P0[2])

        s = (rhos - rho0)/self.dr
        t = (phis - phi0)/self.dphi
        u = (zetas - z0)/self.dz

        return s,t,u 


class FFDCylindrical:
    def __init__(self, gmsh_model, l, m , n, dim, tag=-1, includeBoundary=False, parametric=True):
        """This class generates a cylindrical FFD Lattice 

        Args:
            gmsh_model (gmsh.model): gmsh model (after meshing) 
            l (int): number of points in the x direction
            m (int): number of points in the y direction
            n (int): number of points in the z direction
            dim (int): dimension of gmsh entity - 2 or 3
            tag (int, optional): tag of the entity, -1 returns all tags. Defaults to -1.
        """
        self.l = l
        self.m = m
        self.n = n

        self.Px = np.zeros((l,m,n))
        self.Py = np.zeros((l,m,n))
        self.Pz = np.zeros((l,m,n))

        nodes, coords, param = gmsh_model.mesh.getNodes(dim, tag, includeBoundary, parametric)

        self.base_coords = coords

        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        self.dx = max(xs)-min(xs)
        self.dy = max(ys)-min(ys)
        self.dz = max(zs)-min(zs)

        rhos, phis, zetas = cart2cyl(xs, ys, zs)
        if dim==2:
            self.dr = max(rhos)
        elif dim==3:
            self.dr = max(rhos)
        print(self.dr)

        self.dphi = 2*np.pi/self.m
        self.dz = max(zs)-min(zs)

        for i in range(l):
            for j in range(m):
                for k in range(n):
                    rho = 0  + self.dr * i / (l - 1)
                    phi = self.dphi* j
                    self.Px[i, j, k] = rho * np.cos(phi)
                    self.Py[i, j, k] = rho * np.sin(phi)
                    self.Pz[i, j, k] = min(zs)  + self.dz * k / (n - 1)

        self.P0 = np.array([0, 0, 0])
        
        self.Pr, self.Pphi, self.Pz = cart2cyl(self.Px, self.Py, self.Pz)

    def write_ffd_points(self, name):
        """writes the FFD points as a vtu file (readable using ParaView).
        """
        
        x_ffd = self.Px.flatten()
        y_ffd = self.Py.flatten()
        z_ffd = self.Pz.flatten()
        pointsToVTK(name, x_ffd, y_ffd, z_ffd)
        print("FFD points are saved.")
    
    # def calcSTU(self, coords):
    #     """
    #     Calc STU (parametric) coordinates
    #     """
    #     xs = coords[0::3]
    #     ys = coords[1::3]
    #     zs = coords[2::3]

    #     rhos, phis, zetas = cart2cyl(xs, ys, zs)
    #     print(self.dr)

    #     s = rhos/self.dr
    #     t = phis/(2*np.pi)
    #     u = zetas/self.dz
    #     print(s)

    #     return s,t,u 

    def calcSTU(self, coords):
        """
        Calc STU (parametric) coordinates
        """
        xs = coords[0::3]
        ys = coords[1::3]
        zs = coords[2::3]

        s = (xs - self.P0[0])/self.dx
        t = (ys - self.P0[1])/self.dy
        u = (zs - self.P0[2])/self.dz

        return s,t,u

def get_meshdata(gmsh_model):
    """ Calculates the current mesh data which inputs the deformation function.

    Args:
        gmsh_model (_type_): gmsh.model

    Returns:
        dictionary: mesh data
    """
    mesh_data = {}
    for e in gmsh_model.getEntities():
        mesh_data[e] = (gmsh_model.getBoundary([e]),
                gmsh_model.mesh.getNodes(e[0], e[1]),
                gmsh_model.mesh.getElements(e[0], e[1]))
    return mesh_data

def calcSTU(coords, P0, dx, dy, dz):
    """
    Calc STU (parametric) coordinates
    """
    xs = coords[0::3]
    ys = coords[1::3]
    zs = coords[2::3]

    s = (xs - P0[0])/dx
    t = (ys - P0[1])/dy
    u = (zs - P0[2])/dz

    return s,t,u

def deform_mesh(gmsh_model, mesh_data, FFDVolume):

    for e in mesh_data:
    
        if len(mesh_data[e][1][1])==3:
            old_coord = mesh_data[e][1][1]
            s,t,u = calcSTU(old_coord,FFDVolume.P0, FFDVolume.dx, FFDVolume.dy, FFDVolume.dz)
            Xdef = np.zeros((1,3))

        else:
            old_coords = mesh_data[e][1][1]
            s,t,u = calcSTU(old_coords,FFDVolume.P0, FFDVolume.dx, FFDVolume.dy, FFDVolume.dz)
            old_coords = old_coords.reshape(-1,3)
            Xdef = np.zeros((len(old_coords),3))

        for point, param_s in enumerate(s):
            for i in range(FFDVolume.l):
                for j in range(FFDVolume.m):
                    for k in range(FFDVolume.n):
                        Xdef[point] +=  comb(FFDVolume.l-1,i)*np.power(1-s[point], FFDVolume.l-1-i)*np.power(s[point],i) * \
                                        comb(FFDVolume.m-1,j)*np.power(1-t[point], FFDVolume.m-1-j)*np.power(t[point],j) * \
                                        comb(FFDVolume.n-1,k)*np.power(1-u[point], FFDVolume.n-1-k)*np.power(u[point],k) * \
                                        np.asarray([FFDVolume.Px[i,j,k], FFDVolume.Py[i,j,k], FFDVolume.Pz[i,j,k]])

    new_coord = Xdef.flatten()
        
    gmsh_model.addDiscreteEntity(e[0], e[1], [b[1] for b in mesh_data[e][0]])
    gmsh_model.mesh.addNodes(e[0], e[1], mesh_data[e][1][0], new_coord)
    gmsh_model.mesh.addElements(e[0], e[1], mesh_data[e][2][0], mesh_data[e][2][1], mesh_data[e][2][2])

    return gmsh_model