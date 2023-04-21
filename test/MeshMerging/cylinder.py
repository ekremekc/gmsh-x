import gmsh
import sys

gmsh.initialize()

L = 1
W = 1
H = 1

R = 0.25
gdim = 3

cube = gmsh.model.occ.addCylinder(L,W/2,H/2, L, 0, 0, R, tag=1)
gmsh.model.occ.synchronize()

surfaces = gmsh.model.getEntities(dim=2)
print(surfaces)

volumes = gmsh.model.getEntities(dim=gdim)
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], 1)


gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
gmsh.model.mesh.generate(gdim)

gmsh.write("{}.msh".format("cylinder"))

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
    
gmsh.finalize()
