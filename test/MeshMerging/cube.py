import gmsh
import sys

gmsh.initialize()

L = 1
W = 1
H = 1
gdim = 3

cube = gmsh.model.occ.addBox(0,0,0, L, W, H, tag=1)
gmsh.model.occ.synchronize()

volumes = gmsh.model.getEntities(dim=gdim)
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], 1)


gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
#gmsh.model.mesh.generate(gdim)

gmsh.write("{}.msh".format("cube"))

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
    
gmsh.finalize()
