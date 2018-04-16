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

# class Front(dolfin.SubDomain):
#     def inside(self, x, on_boundary):
#         return dolfin.near(x[2], 0.0)
# 
# class Back(dolfin.SubDomain):
#     def inside(self, x, on_boundary):
#         return dolfin.near(x[2], 1.0)

# INTERNAL BOUNDARY CONDITION
# ELECTRODE
class Obstacle(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        # return (dolfin.between(x[0], (0.4, 0.6)) and dolfin.between(x[1], (0.1, 0.6)) and dolfin.between(x[2], (0.4, 0.6)))
        return (dolfin.between(x[0], (0.4, 0.6)) and dolfin.between(x[1], (0.4, 0.6)) )

# INTERNAL BOUNDARY CONDITION
#DIELECTRIC
class Obstacle2(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        # return (dolfin.between(x[0], (0.3, 0.7)) and dolfin.between(x[1], (0.7, 0.9)) and dolfin.between(x[2], (0.3, 0.7)))
        return (dolfin.between(x[1], (0.8, 1.0)) )

# Initialize sub-domain instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()
# front = Front()
# back = Back()
obstacle = Obstacle()
obstacle2 = Obstacle2()



# Define mesh

domain = '''
//mesh_resolution: small->fine, large->coarse
mesh_resolution = 0.01;
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
Point(5) = {0.01, 0.2, 0, mesh_resolution};
Point(6) = {0.99, 0.2, 0, mesh_resolution}; 
Line(7) = {5, 6};
Physical Line(1) = {7};
Line{7} In Surface{6};

// 2
Point(7) = {0.01, 0.8, 0, mesh_resolution};
Point(8) = {0.99, 0.8, 0, mesh_resolution};
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
Physical Line(5) = {12};
Line{12} In Surface{6};

Physical Surface(1) = {6};
'''

with open('/home/nathan/git/PoisSolvers/FEM/data/domain.geo', 'w') as f: f.write(domain)
subprocess.call(['gmsh -2 /home/nathan/git/PoisSolvers/FEM/data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert /home/nathan/git/PoisSolvers/FEM/data/domain.msh domain.xml'], shell=True)
mesh = dolfin.Mesh('domain.xml')

# domain = mshr.Box(dolfin.Point(0,0,0), dolfin.Point(1,1,1))
# mesh = mshr.generate_mesh(domain, 50)

print "Mesh initialized"
# Initialize mesh function for interior domains
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
obstacle2.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
# front.mark(boundaries, 5)
# back.mark(boundaries, 6)
obstacle.mark(boundaries, 5)
print "Boundaries marked"

# Define input data
a0 = dolfin.Constant(1.0)
a1 = dolfin.Constant(100.0)
# g_L = Expression("- 10*exp(- pow(x[1] - 0.5, 2))", degree=2)
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
f = dolfin.Constant(0.05)

print "Inputs finished"

# Define function space and basis functions
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

print "Function, space, basis defined"
# Define Dirichlet boundary conditions at top, bottom, front, back
bcs = [dolfin.DirichletBC(V, 0.0, boundaries, 1),
       dolfin.DirichletBC(V, 0.0, boundaries, 3),
       # dolfin.DirichletBC(V, 0.0, boundaries, 5),
       # dolfin.DirichletBC(V, 0.0, boundaries, 6),
       dolfin.DirichletBC(V, 20, boundaries, 5)]

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
dolfin.solve(a == L, u, bcs, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'ilu'})
# solve(a == L, u, bcs, solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})
print "Solve finished"
# Evaluate integral of normal gradient over top boundary

# Create boundary mesh consisting of the outside
bmesh = dolfin.BoundaryMesh(mesh, "exterior")

n = dolfin.FacetNormal(mesh)
m1 = dolfin.dot(dolfin.grad(u), n)*dolfin.ds(2)
v1 = dolfin.assemble(m1)
print "\int grad(u) * n ds(2) = ", v1

# Evaluate integral of u over the obstacle
m2 = u*dolfin.dx(1)
v2 = dolfin.assemble(m2)
print "\int u dx(1) = ", v2

print "Print solution to file."
file_string = "/home/nathan/git/PoisSolvers/FEM/data/Potential.pvd"
dolfin.File(file_string) << u

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
