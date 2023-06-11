import gmsh
import sys
import numpy as np

gmsh.initialize()

r_outer = 0.25
z_inlet = 0.
z_outlet = 0.5
lc = 0.01

deg = 2

N = 40
angle_jump = 2*np.pi/N

N_wall = 3
z_wall = np.linspace(z_inlet,z_outlet,N_wall)

# inlet
angle = 0

for i in range(N+1): 
    for z_point in z_wall:  
        if z_point==0.25:
            gmsh.model.occ.addPoint(r_outer*np.cos(angle), r_outer*np.sin(angle)+0.4, z_point,lc)
        else:
            gmsh.model.occ.addPoint(r_outer*np.cos(angle), r_outer*np.sin(angle), z_point,lc)
    angle += angle_jump

gmsh.model.occ.addBSplineSurface(range(1,(N+1)*N_wall+1), N_wall)

gmsh.model.occ.removeAllDuplicates()

# inlet surface
tolerance = 0.01
inlet_points = gmsh.model.occ.getEntitiesInBoundingBox(-r_outer-tolerance,-r_outer-tolerance,z_inlet-tolerance,
                                                        r_outer+tolerance,r_outer+tolerance,z_inlet+tolerance,
                                                        0)
inlet_points_tags = [point[1] for point in inlet_points]
# print(inlet_points_tags)

C_inlet = gmsh.model.occ.addBSpline(inlet_points_tags, degree=deg)
W_inlet = gmsh.model.occ.addWire([C_inlet])
gmsh.model.occ.addSurfaceFilling(W_inlet)


# outlet surface
outlet_points = gmsh.model.occ.getEntitiesInBoundingBox(-r_outer-tolerance,-r_outer-tolerance,z_outlet-tolerance,
                                                        r_outer+tolerance,r_outer+tolerance,z_outlet+tolerance,
                                                        0)
outlet_points_tags = [point[1] for point in outlet_points]
# print(outlet_points_tags)

C_outlet = gmsh.model.occ.addBSpline(outlet_points_tags, degree=deg)
W_outlet = gmsh.model.occ.addWire([C_outlet])
gmsh.model.occ.addSurfaceFilling(W_outlet)

# Volume of the pipe
gmsh.model.occ.removeAllDuplicates()

gmsh.model.occ.synchronize()

# surfaces = gmsh.model.occ.getEntities(2)
# surface_tags = [surface[1] for surface in surfaces]
# print(surface_tags)
gmsh.model.occ.remove([(2, 2)])
gmsh.model.occ.remove([(2, 4)])
gmsh.model.occ.addSurfaceLoop([1,3,5], 1)
gmsh.model.occ.addVolume([1], 1)

gmsh.model.occ.synchronize()

# Add physical tags

gmsh.model.addPhysicalGroup(2,[1],1) # Lateral
gmsh.model.addPhysicalGroup(2,[3],2) # Inlet
gmsh.model.addPhysicalGroup(2,[5],3) # Outlet
gmsh.model.addPhysicalGroup(1,[1],1) # Pipe

gmsh.option.setNumber('Mesh.MeshSizeMax', 0.02)

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()