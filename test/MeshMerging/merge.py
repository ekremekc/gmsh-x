import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.merge("cube.msh")
gmsh.merge("cylinder.msh")
gmsh.model.occ.synchronize()

surfaces = gmsh.model.getEntities(dim=2)
print("Surfaces:", surfaces)

volumes = gmsh.model.getEntities(dim=3)
print("Volumes:", volumes)
gmsh.write("{}.msh".format("merged"))

gmsh.fltk.run()
gmsh.finalize()
