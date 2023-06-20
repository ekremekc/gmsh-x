import sys
import numpy as np
import gmsh

geom_dir = "GeomDir/"
mesh_dir = "MeshDir/"

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("cubes")
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

gmsh.model.occ.importShapes(geom_dir+'cube1.step')
gmsh.model.occ.synchronize()
gmsh.model.occ.importShapes(geom_dir+'cube2.step')
gmsh.model.occ.synchronize()
gmsh.model.occ.importShapes(geom_dir+'cube3.step')
gmsh.model.occ.synchronize()

gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

gmsh.model.addPhysicalGroup(3, [1], tag=1) 
gmsh.model.addPhysicalGroup(3, [2], tag=2) 
gmsh.model.addPhysicalGroup(3, [3], tag=3) 


surfaces = gmsh.model.occ.getEntities(dim=2)

for surface in surfaces:

    gmsh.model.addPhysicalGroup(surface[0], [surface[1]])

gmsh.model.occ.synchronize()
lc = 0.01
gmsh.model.mesh.field.add("Constant", 1)
gmsh.model.mesh.field.setNumbers(1, "VolumesList", [2])
gmsh.model.mesh.field.setNumber(1, "VIn", lc / 8)
gmsh.model.mesh.field.setNumber(1, "VOut", lc)


gmsh.model.mesh.field.add("Constant", 2)
gmsh.model.mesh.field.setNumbers(2, "VolumesList", [3])
gmsh.model.mesh.field.setNumber(2, "VIn", lc / 20)
gmsh.model.mesh.field.setNumber(2, "VOut", lc)

gmsh.model.mesh.field.add("Min", 3)
gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1,2])


gmsh.model.mesh.field.setAsBackgroundMesh(3)

gmsh.model.occ.synchronize()

#gmsh.option.setNumber("Mesh.MeshSizeMin", 0.005)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.01)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)#10
gmsh.option.setNumber("Mesh.RandomFactor", 1e-11)
gmsh.option.setNumber("Mesh.RandomFactor3D", 1e-13)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
gmsh.model.mesh.generate(3)

gmsh.write("{}.msh".format(mesh_dir+"cubes"))

#if '-nopopup' not in sys.argv:
#    gmsh.fltk.run()

gmsh.finalize()

from helmholtz_x.dolfinx_utils import write_xdmf_mesh
write_xdmf_mesh(mesh_dir+"cubes",dimension=3)
    

