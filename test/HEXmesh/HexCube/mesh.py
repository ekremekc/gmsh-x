import sys
import numpy as np
import gmsh

mesh_dir = "MeshDir/"
geom_dir = "GeomDir/"
geom_name = "cube"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

gmsh.model.add("cube")
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

gmsh.model.occ.importShapes(geom_dir+geom_name+".step")
gmsh.model.occ.synchronize()


gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(3, [1], tag=1) 


surfaces = gmsh.model.occ.getEntities(dim=2)

for surface in surfaces:

    gmsh.model.addPhysicalGroup(surface[0], [surface[1]])

gmsh.model.occ.synchronize()


gmsh.model.occ.synchronize()

transfinite = True
transfiniteAuto = False

if transfinite:
    NN = 20
    for c in gmsh.model.getEntities(1):
        print(c)
        gmsh.model.mesh.setTransfiniteCurve(c[1], NN)
    for s in gmsh.model.getEntities(2):
        print("Surface: ", s)
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(s[0], s[1])
        gmsh.model.mesh.setSmoothing(s[0], s[1], 10)
    gmsh.model.mesh.setTransfiniteVolume(1)
    
elif transfiniteAuto:
    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.1)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.1)
    # setTransfiniteAutomatic() uses the sizing constraints to set the number
    # of points
    gmsh.model.mesh.setTransfiniteAutomatic()

else:
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)


gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)#10
gmsh.option.setNumber("Mesh.RandomFactor", 1e-11)
gmsh.option.setNumber("Mesh.RandomFactor3D", 1e-13)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
gmsh.model.mesh.generate(3)

gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.write("{}.msh".format(mesh_dir+geom_name))

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

from helmholtz_x.dolfinx_utils import write_xdmf_hexmesh, XDMFReader
write_xdmf_hexmesh(mesh_dir+geom_name,dimension=3)

geom = XDMFReader(mesh_dir+geom_name)
geom.getInfo()

