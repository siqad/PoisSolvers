# from dolfin import *
import dolfin
import matplotlib.pyplot as plt
import mshr
from paraview import simple as pvs

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

class Front(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], 0.0)

class Back(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], 1.0)


# INTERNAL BOUNDARY CONDITION
class Obstacle(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[0], (0.4, 0.6)) and dolfin.between(x[1], (0.1, 0.9)) and dolfin.between(x[2], (0.4, 0.6)))

# Initialize sub-domain instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()
front = Front()
back = Back()
obstacle = Obstacle()

# Define mesh
domain = mshr.Box(dolfin.Point(0,0,0), dolfin.Point(1,1,1))
mesh = mshr.generate_mesh(domain, 25)
# mesh = dolfin.UnitCubeMesh(10, 10, 10)
# mesh = UnitSquareMesh(64, 64)
print "Mesh initialized"
# Initialize mesh function for interior domains
# domains = CellFunction("size_t", mesh)
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
obstacle.mark(domains, 1)

# Initialize mesh function for boundary domains
# boundaries = FacetFunction("size_t", mesh)
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
front.mark(boundaries, 5)
back.mark(boundaries, 6)

print "Boundaries marked"

# Define input data
a0 = dolfin.Constant(1.0)
a1 = dolfin.Constant(100.0)
# g_L = Expression("- 10*exp(- pow(x[1] - 0.5, 2))", degree=2)
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("1.0")
f = dolfin.Constant(0.05)

print "Inputs finished"

# Define function space and basis functions
# V = FunctionSpace(mesh, "CG", 2)
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

print "Function, space, basis defined"
# Define Dirichlet boundary conditions at top and bottom boundaries
bcs = [dolfin.DirichletBC(V, 0.0, boundaries, 1),
       dolfin.DirichletBC(V, 0.0, boundaries, 3),
       dolfin.DirichletBC(V, 0.0, boundaries, 5),
       dolfin.DirichletBC(V, 0.0, boundaries, 6)]

print "Dirichlet boundaries defined"
# Define new measures associated with the interior domains and
# exterior boundaries

# dx = Measure("dx")[domains]
dx = dolfin.dx(subdomain_data=domains)
# ds = Measure("ds")[boundaries]
ds = dolfin.ds(subdomain_data=boundaries)
print "Measures defined"

# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    + dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(2) \
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
# # Evaluate integral of normal gradient over top boundary
#
# # Create boundary mesh consisting of the 6 sides of the cube
# bmesh = dolfin.BoundaryMesh(mesh, "exterior")
#
# # Create SubMesh for side at z=0
# # This will be a UnitSquareMesh with topology dimension 2 in 3 space dimensions
#
# cc = dolfin.MeshFunction("size_t", bmesh, bmesh.topology().dim())
# # cc = CellFunction('size_t', bmesh, 0)
# xyplane = dolfin.AutoSubDomain(lambda x: x[2] < 1e-8)
# xyplane.mark(cc, 1)
# submesh = dolfin.SubMesh(bmesh, cc, 1)
#
# # Move slice/submesh to z=0.5
# x = submesh.coordinates()
# x[:, 2] += 0.5
#
# # plt.figure()
# # plot(submesh, title="submesh")
#
# # Create a FunctionSpace on the submesh
# Vs = dolfin.FunctionSpace(submesh, "CG", 1)
# u.set_allow_extrapolation(True)
#
# # interpolate_nonmatching_mesh required in parallel,
# # interpolate works in series
# # us = interpolate_nonmatching_mesh(u, Vs)
# us = dolfin.interpolate(u, Vs)
#
# plt.figure()
# dolfin.plot(us, title="slice")
#
# plt.figure()
# dolfin.plot(dolfin.grad(us), title="Projected grad(us)")

n = dolfin.FacetNormal(mesh)
m1 = dolfin.dot(dolfin.grad(u), n)*dolfin.ds(2)
v1 = dolfin.assemble(m1)
print "\int grad(u) * n ds(2) = ", v1

# Evaluate integral of u over the obstacle
m2 = u*dolfin.dx(1)
v2 = dolfin.assemble(m2)
print "\int u dx(1) = ", v2

# Plot solution and gradient
# plt.figure()
# dolfin.plot(u, title="u")
# plt.figure()
# dolfin.plot(dolfin.grad(u), title="Projected grad(u)")
# plt.figure()
# dolfin.plot(mesh, title="Mesh")


# plt.show()
# interactive()

print "Print solution to file."
dolfin.File("u.pvd") << u

##PLOTTING
data_3d = pvs.PVDReader(FileName="u.pvd")
#create the slice
slice_2d = pvs.Slice(Input=data_3d, SliceType="Plane")
slice_2d.SliceType.Origin = [0, 0, 0.5]
slice_2d.SliceType.Normal = [0, 0, 1]
pvs.Show(data_3d)
prop = pvs.GetDisplayProperties(data_3d)
prop.Representation = "Wireframe"
print prop.GetProperty("Representation")
pvs.Show(slice_2d)
pvs.Interact(view=None)
pvs.Render()
