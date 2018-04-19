# from dolfin import *
import dolfin
import matplotlib.pyplot as plt
import mshr
from paraview import simple as pvs
import subprocess

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], 0.0)

class Right(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], 1.0)

class Bottom(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], 0.0)

class Top(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], 1.0)

# # INTERNAL BOUNDARY CONDITION
# ELECTRODE
class Obstacle(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return not on_boundary and (dolfin.between(x[0], (0.4, 0.6)) and dolfin.between(x[1], (0.4,0.6))) 
        # return (dolfin.near(x[0], 0.4) or dolfin.near(x[0], 0.6)) and (dolfin.near(x[1], 0.4) or dolfin.near(x[1], 0.6)) 

# INTERNAL BOUNDARY CONDITION
#DIELECTRIC
class Obstacle2(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[1], (0.8, 1.0)))

# Initialize sub-domain instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()
obstacle = Obstacle()
obstacle2 = Obstacle2()

# Define mesh
domain = '''
//mesh_resolution: small->fine, large->coarse
mesh_resolution = 0.05;
Point(1) = {0, 0, 0, mesh_resolution};
Point(2) = {1, 0, 0, mesh_resolution};
Point(3) = {1, 1, 0, mesh_resolution};
Point(4) = {0, 1, 0, mesh_resolution};
// 2d domain
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Cracks are locations that the mesh is forced to align to

// 1
//Point(5) = {0.01, 0.2, 0, mesh_resolution};
//Point(6) = {0.99, 0.2, 0, mesh_resolution}; 
//Line(7) = {5, 6};
//Physical Line(1) = {7};
//Line{7} In Surface{6};

// 2
Point(7) = {0.1, 0.8, 0, mesh_resolution};
Point(8) = {0.9, 0.8, 0, mesh_resolution};
Line(8) = {7, 8};
Physical Line(2) = {8};
Line{8} In Surface{6};


Point(9) = {0.4, 0.4, 0, mesh_resolution};
Point(10) = {0.4, 0.6, 0, mesh_resolution};
Point(11) = {0.6, 0.4, 0, mesh_resolution};
Point(12) = {0.6, 0.6, 0, mesh_resolution};

Line(9) = {9, 10};
Physical Line(3) = {9};
Line{9} In Surface{6};
Line(10) = {10, 12};
Physical Line(4) = {10};
Line{10} In Surface{6};
Line(11) = {12, 11};
Physical Line(5) = {11};
Line{11} In Surface{6};
Line(12) = {11, 9};
Physical Line(6) = {12};
Line{12} In Surface{6};

Physical Surface(1) = {6};
'''

with open('../data/domain.geo', 'w') as f: f.write(domain)
subprocess.call(['gmsh -2 ../data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
mesh = dolfin.Mesh('domain.xml')

# domain = mshr.Rectangle(dolfin.Point(0,0), dolfin.Point(1,1))
# mesh = mshr.generate_mesh(domain, 20)

# mesh = dolfin.RectangleMesh(dolfin.Point(0.0,0.0),dolfin.Point(1.0,1.0),20,20)
print "Mesh initialized"
# Initialize mesh function for interior domains
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
# obstacle2.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
obstacle.mark(boundaries, 5)
print "Boundaries marked"

# Define input data
a0 = dolfin.Constant(1.0)
a1 = dolfin.Constant(100.0)
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
f = dolfin.Constant(0.10)

print "Inputs finished"

# Define function space and basis functions
V = dolfin.FunctionSpace(mesh, "CG", 2)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

print "Function, space, basis defined"
# Define Dirichlet boundary conditions at top, bottom, front, back
bcs = [dolfin.DirichletBC(V, 0.0, boundaries, 1),
       dolfin.DirichletBC(V, 0.0, boundaries, 3),
       dolfin.DirichletBC(V, 20.0, boundaries, 5)]

print "Dirichlet boundaries defined"
# Define new measures associated with the interior domains and
# exterior boundaries

dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)
print "Measures defined"

# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    # + dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(2) \
    - g_L*v*ds(2) - g_R*v*ds(4) \
    - f*v*dx(0) - f*v*dx(1) ) #- f*v*dolfin.dx(2) )
print "Variational Form defined"
# Separate left and right hand sides of equation
a, L = dolfin.lhs(F), dolfin.rhs(F)
print "Separated LHS and RHS"
# Solve problem
u = dolfin.Function(V)
M = u*dx()
tol = 1.0e-8

dolfin.parameters["refinement_algorithm"] = "plaza_with_parent_facets"
problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
solver.solve()
# solver = dolfin.AdaptiveLinearVariationalSolver(problem, M)
# solver.solve(tol)
# solver.summary()

# dolfin.solve(a == L, u, bcs, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'ilu'})
# solve(a == L, u, bcs, solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})
print "Solve finished"

print "Print solution to file."
file_string = "../data/Potential.pvd"
dolfin.File(file_string) << u

plt.figure()
dolfin.plot(u, title="u")
# plt.figure()
# dolfin.plot(u.root_node(), title="Solution on initial mesh")
# plt.figure()
# dolfin.plot(u.leaf_node(), title="Solution on final mesh")

plt.show()
##PLOTTING
data_3d = pvs.PVDReader(FileName=file_string)
#create the slice
# slice_2d = pvs.Slice(Input=data_3d, SliceType="Plane")
# slice_2d.SliceType.Origin = [0, 0, 0.5]
# slice_2d.SliceType.Normal = [0, 0, 1]
pvs.Show(data_3d)
prop = pvs.GetDisplayProperties(data_3d)
prop.Representation = "Wireframe"
print prop.GetProperty("Representation")
# pvs.Show(slice_2d)
pvs.Interact(view=None)
pvs.Render()
