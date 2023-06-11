from gmsh_x.nurbs import nurbs_annulus, nurbs_cylinder, generate_radius_dict
from helmholtz_x.dolfinx_utils import dict_loader
import sys
import gmsh
import numpy as np

step = 0.1
# Mesh generation
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Plenum geometry
r_outer = 0.21
r_inner = 0.14

z_start = 0.
length = 0.4
z_end = z_start + length

lc = length/10

N_u = 9   # number of control points in the circumferential direction
N_cpt = 5 # number of control points in the +z direction 

r_outer = generate_radius_dict(r_outer, N_u, N_cpt)

# # Inner Plenum
r_inner = generate_radius_dict(r_inner, N_u, N_cpt)

annulus = nurbs_annulus(z_start, z_end, r_outer, r_inner, N_cpt, lc)

# gmsh.option.setNumber('Mesh.MeshSizeMin', 0.01) 
# gmsh.option.setNumber('Mesh.MeshSizeMax', 0.03) 

gmsh.model.mesh.generate(3)
# gmsh.model.mesh.refine()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()