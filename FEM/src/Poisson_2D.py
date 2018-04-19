# from dolfin import *
import dolfin
import matplotlib.pyplot as plt
import mshr
from paraview import simple as pvs
import subprocess
import mesh_writer

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], 0.0)

class Right(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], 10.0)

class Bottom(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], 0.0)

class Top(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], 10.0)

# INTERNAL BOUNDARY CONDITION
# ELECTRODE
class Electrode(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[0], (4.0, 6.0)) and dolfin.between(x[1], (4.0, 6.0)) )

# INTERNAL BOUNDARY CONDITION
#DIELECTRIC
class Dielectric(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[1], (8.0, 10.0)) )

# Initialize sub-domain instances
left = Left()
top = Top()
right = Right()
bottom = Bottom()
electrode = Electrode()
dielectric = Dielectric()

mw = mesh_writer.MeshWriter()
mw.resolution = 1.0
mw.addOuterBound([0.0,0.0], [10.0,10.0], 1)
mw.addCrack([0.01,8.0],[9.99,8.0],0.1)
mw.addCrack([0.01,7.0],[9.99,7.0],0.5)
mw.addCrackBox([4.0,4.0],[6.0,6.0],0.1)
# print mw.file_string

with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
subprocess.call(['gmsh -2 ../data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
mesh = dolfin.Mesh('domain.xml')

print "Mesh initialized"
# Initialize mesh function for interior domains
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
dielectric.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
# front.mark(boundaries, 5)
# back.mark(boundaries, 6)
electrode.mark(boundaries, 5)
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
print "Solve finished"

print "Print solution to file."
file_string = "../data/Potential.pvd"
dolfin.File(file_string) << u

plt.figure()
dolfin.plot(u, title="u")
dolfin.plot(mesh)
plt.show()
##PLOTTING
data_3d = pvs.PVDReader(FileName=file_string)
#create the slice
# slice_2d = pvs.Slice(Input=data_3d, SliceType="Plane")
# slice_2d.SliceType.Origin = [0, 0, 0.5]
# slice_2d.SliceType.Normal = [0, 0, 1]
# pvs.Show(data_3d)
# prop = pvs.GetDisplayProperties(data_3d)
# prop.Representation = "Wireframe"
# print prop.GetProperty("Representation")
# # pvs.Show(slice_2d)
# pvs.Interact(view=None)
# pvs.Render()
