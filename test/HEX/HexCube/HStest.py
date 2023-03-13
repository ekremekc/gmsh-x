from helmholtz_x.eigensolvers_x import eps_solver
from helmholtz_x.passive_flame_x import PassiveFlame
from helmholtz_x.eigenvectors_x import normalize_eigenvector
from helmholtz_x.dolfinx_utils import xdmf_writer, XDMFReader
from helmholtz_x.solver_utils import start_time, execution_time
from helmholtz_x.parameters_utils import c_uniform
start = start_time()
import numpy as np

# approximation space polynomial degree
degree = 1

# number of elements in each direction of mesh
cube = XDMFReader("MeshDir/cube")
mesh, subdomains, facet_tags = cube.getAll()
cube.getInfo()

# Define the boundary conditions
boundary_conditions = {}

# Define Speed of sound

sos = 340
c = c_uniform(mesh, sos)

# Introduce Passive Flame Matrices

matrices = PassiveFlame(mesh, subdomains, facet_tags, boundary_conditions, c, degree=degree)

matrices.assemble_A()
matrices.assemble_C()

target = 200 * 2 * np.pi 
E = eps_solver(matrices.A, matrices.C, target**2, nev=2, print_results= True)

omega, uh = normalize_eigenvector(mesh, E, 1, absolute=True, degree=degree, which='right')

xdmf_writer("Results/P", mesh, uh)

execution_time(start)
