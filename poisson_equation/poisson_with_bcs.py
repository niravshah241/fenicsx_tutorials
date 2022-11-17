import gmsh
import os
import numpy as np
import matplotlib.pyplot as plt

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx.mesh import *
from dolfinx.fem import *
from dolfinx.io import *
from ufl import *

import ufl
import dolfinx

gmsh.initialize('',False) # Avoid reading config file such as .gmshrc as it may have old data

gdim = 2 # Geometric dimensions of the domain

h1 = 0.2 # Very coarse
h2 = 0.1 # Coarse enough
h3 = 0.05 # Fine enough
h4 = 0.02 # Very fine

H = 10. # Height of square
W = 10. # Width of square
r = 5. # radius of circle

assert r < H
assert r < W

# Points
A = gmsh.model.occ.addPoint(0., 0., 0., h1)
B = gmsh.model.occ.addPoint(0., H, 0., h1)
C = gmsh.model.occ.addPoint(W, H, 0., h1)
D = gmsh.model.occ.addPoint(W, 0., 0., h1)

# Lines
AB = gmsh.model.occ.addLine(A,B)
BC = gmsh.model.occ.addLine(B,C)
CD = gmsh.model.occ.addLine(C,D)
DA = gmsh.model.occ.addLine(D,A)

# Surfaces and final geometry
ABCDA = gmsh.model.occ.addCurveLoop([AB,BC,CD,DA])
rectangle = gmsh.model.occ.addPlaneSurface([ABCDA],1) # rectangle
c_r = gmsh.model.occ.addCircle(0., H, 0., r) # Creates Lines "equivalent" of circumference
c_r_curve = gmsh.model.occ.addCurveLoop([c_r],2) # Creates curve of circumference
circle = gmsh.model.occ.addPlaneSurface([c_r_curve],2) # Circle (plane surface type)
#circle = gmsh.model.occ.addDisk(0, H, 0, r, r) # NOTE Adding disk is equivalent to creating circle from above three lines
domain_geometry = gmsh.model.occ.cut([(gdim,rectangle)],[(gdim,circle)],3) # 3rd argument is marker of the surface domain_geometry

# Synchronize
gmsh.model.occ.synchronize()

# Create mesh
gmsh.option.setNumber("Mesh.Algorithm", 8) # 8=Frontal-Delaunay for Quads (See section 7.4,  https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2) # 2=simple full-quad (See section 7.4,  https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options)
gmsh.option.setNumber("Mesh.RecombineAll", 1) # Apply recombination algorithm to all surfaces, ignoring per-surface spec
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) # Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3: barycentric)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1) # Minimum characteristic element size
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.3) # Maximum characteristic element size
gmsh.model.mesh.generate(gdim) # Mesh generation
gmsh.model.mesh.setOrder(1) # Mesh order
gmsh.model.mesh.optimize("Netgen") # Mesh optimisation or improving quality of mesh

# Extract edges and surfaces to add physical groups
surfaces = gmsh.model.getEntities(dim=gdim)
edges = gmsh.model.getBoundary(surfaces) # Gives 'list' of boundaries in the form [(gdim-1),marker] with length = number of boundaries
# edges = gmsh.model.getEntities(dim=gdim-1) # Alternate to getting boundaries

for i in range(1,len(surfaces)+1):
    gmsh.model.addPhysicalGroup(gdim,[surfaces[i-1][1]],surfaces[i-1][1])
for i in range(1,len(edges)+1):
    gmsh.model.addPhysicalGroup(gdim-1,[edges[i-1][1]],edges[i-1][1])

# Import mesh in dolfinx
gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)

assert domain.topology.dim == gdim, "Geometric dimensions of mesh and domain are not same"

with XDMFFile(MPI.COMM_WORLD, "mesh_data/mesh.xdmf", "w") as mesh_file_xdmf:
    mesh_file_xdmf.write_mesh(domain)
    mesh_file_xdmf.write_meshtags(facet_markers)

V_cg2 = FiniteElement("CG",domain.ufl_cell(),2)
V = dolfinx.fem.FunctionSpace(domain,V_cg2) #NOTE Only FunctionSpace conflicts with ufl.FunctionSpace and dolfinx.fem.FunctionSpace . Always use fem.FunctionSpace in python. 
V2 = dolfinx.fem.FunctionSpace(domain,("CG",2)) #NOTE To be used for analytical solution during error analysis
u, v = TrialFunction(V), TestFunction(V)

def u_analytical(x):
    return 1 + x[0]**2 + 2*x[1]**2

def u_bc_4(x):
    values = np.zeros((1,x.shape[1]), dtype = PETSc.ScalarType) # NOTE Shape of np.zeros is dim_of_physical_quantity \times number of nodes. PETSc.ScalarType is necessary to increase portability and proper data type in interpolation. 
    values[0] = 1 + x[0]**2 + 2*x[1]**2
    return values

x = ufl.SpatialCoordinate(domain) # TODO find alternative to using SpatialCoordinate
A = PETSc.ScalarType(1.)
B = PETSc.ScalarType(2.)
C = PETSc.ScalarType(1.)
u_ufl = C + A * x[0]**2 + B * x[1]**2

''' NOTE
Defining as constant PETSc.ScalarType coefficients is advantageous for repeated compilation of expression. The coefficients are like pointers which are evaluated from a memory location instead of read as fixed numbers.  
'''

class exact_solution():
    '''
    class with __call__ method. Define boundary condition (= exact solution in this case) in __call__ method.
    '''
    def __init__(self):
        print("Initialised the class")
    def __call__(self,x):
        return 1 + x[0]**2 + 2*x[1]**2
u_bc_5 = exact_solution() 

#NOTE Function interpolate in fenicsx takes callable objects such as lambda function, python functions, class with __call__ method returning boundary value. This interpolated function can then be used as an input to dolfinx.fem.dirichletbc. The function interpolation also takes non callable argument as dolfinx.fem.Expression or another dolfinx.fem.Function (even with different meshes (important for multigrid)).
# TODO for later: Why use function defined over domain for BCs. Why not Collapse the function on the boundary and use that in BC? 
uD_1 = Function(V)
uD_1.interpolate(lambda x: 1 + x[0]**2 + 2*x[1]**2) 
uD_2 = Function(V)
uD_2.interpolate(lambda x: eval(str(u_ufl)))
uD_3 = Function(V)
uD_3.interpolate(dolfinx.fem.Expression(u_ufl,V.element.interpolation_points())) #NOTE ok but inefficient/expensive but good for coupling
uD_4 = Function(V)
uD_4.interpolate(u_bc_4)
uD_5 = Function(V)
uD_5.interpolate(u_bc_5)

dofs_bc_1 = locate_dofs_topological(V, gdim-1, facet_markers.find(1))
dofs_bc_2 = locate_dofs_topological(V, gdim-1, facet_markers.find(2))
dofs_bc_3 = locate_dofs_topological(V, gdim-1, facet_markers.find(3))
dofs_bc_4 = locate_dofs_topological(V, gdim-1, facet_markers.find(4))
dofs_bc_5 = locate_dofs_topological(V, gdim-1, facet_markers.find(5))

bc_1 = dirichletbc(uD_1,dofs_bc_1)
bc_2 = dirichletbc(uD_2,dofs_bc_2)
bc_3 = dirichletbc(uD_3,dofs_bc_3)
bc_4 = dirichletbc(uD_4,dofs_bc_4)
bc_5 = dirichletbc(uD_5,dofs_bc_5)

a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
f = -div(grad(u_ufl)) # dolfinx.fem.Constant(domain,PETSc.ScalarType(-6.))
L = ufl.inner(f,v) * ufl.dx
problem = petsc.LinearProblem(a, L, bcs=[bc_1,bc_2,bc_3,bc_4,bc_5], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()
uh.name = "computed_solution" # NOTE No space in function names

u_exact = dolfinx.fem.Function(V2,name="exact_solution")
u_exact.interpolate(u_analytical)

H1_error = form(inner(uh - u_exact, uh - u_exact) * dx + inner(grad(uh - u_exact), grad(uh-u_exact)) * dx)
error_local = assemble_scalar(H1_error)
error_H1 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

print("Error in H1 norm: ",error_H1)

with XDMFFile(MPI.COMM_WORLD, "solution_poisson/solution_u.xdmf", "w") as ufile_xdmf:
    uh.x.scatter_forward() # NOTE Communicate the actual value to ghost in other process.
    ufile_xdmf.write_mesh(domain)
    ufile_xdmf.write_function(uh)
    ufile_xdmf.write_function(u_exact)

gmsh.finalize()

'''
NOTE
Ways to impose constant boundary condition (as an example take bc_1 as constant value of 2. on boundary 1):

1. bc_1 = dirichletbc(PETSc.ScalarType(2.),dofs_bc_1,V)
2. bc_1 = dirichletbc(dolfinx.fem.Constant(domain,PETSc.ScalarType(2.)),dofs_bc_1,V)

Notice that third argument FunctionSpace is added.
'''
