# from dolfin import *
import dolfin
import matplotlib.pyplot as plt
import mshr
from paraview import simple as pvs
import subprocess
import mesh_writer


elec_length = 1.0
elec_spacing = 1.5
boundary_x_min = 0.0
boundary_x_max = 10.0
boundary_y_min = 0.0
boundary_y_max = 10.0
mid_x = (boundary_x_max+boundary_x_min)/2.0
mid_y = (boundary_y_max+boundary_y_min)/2.0

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_x_min
        return dolfin.near(x[0], boundary_x_min)

class Right(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_x_max
        return dolfin.near(x[0], boundary_x_max)

class Bottom(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_y_min
        return dolfin.near(x[1], boundary_y_min)

class Top(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_y_max
        return dolfin.near(x[1], boundary_y_max)

# INTERNAL BOUNDARY CONDITION
# ELECTRODE
class Electrode(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x, mid_y
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length, mid_x+0.5*elec_length)) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) or \
           (dolfin.between(x[0], (mid_x-0.5*elec_length-(elec_length+elec_spacing), mid_x+0.5*elec_length-(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) or \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+(elec_length+elec_spacing), mid_x+0.5*elec_length+(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) or \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+2*(elec_length+elec_spacing), mid_x-0.01+2*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) or \
           (dolfin.between(x[0], (mid_x+0.01-2*(elec_length+elec_spacing), mid_x+0.5*elec_length-2*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) \
        )

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
mw.addCrackBox([mid_x-0.5*elec_length,mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length,mid_y+0.5*elec_length],0.1)
mw.addCrackBox([mid_x-0.5*elec_length-(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length-(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addCrackBox([mid_x-0.5*elec_length+(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length+(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)

mw.addCrackBox([mid_x+0.01-2*(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length-2*(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addCrackBox([mid_x-0.5*elec_length+2*(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x-0.01+2*(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)

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
electrode.mark(boundaries, 5)
print "Boundaries marked"

# Define input data
a0 = dolfin.Constant(1.0)
a1 = dolfin.Constant(100.0)
g_T = dolfin.Constant("-1.0")
g_B = dolfin.Constant("-1.0")
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
bcs = [#dolfin.DirichletBC(V, 0.0, boundaries, 1),
       #dolfin.DirichletBC(V, 0.0, boundaries, 3),
       # dolfin.DirichletBC(V, 0.0, boundaries, 5),
       # dolfin.DirichletBC(V, 0.0, boundaries, 6),
       dolfin.DirichletBC(V, 20, boundaries, 5)]

print "Dirichlet boundaries defined"
# Define new measures associated with the interior domains and
# exterior boundaries

dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)
print "Measures defined"

# 1 left, 2 top, 3 right, 4 bottom
# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    # + dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(2) \
    - g_L*v*ds(1) - g_R*v*ds(3) \
    - g_T*v*ds(2) - g_B*v*ds(4) \
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
# dolfin.plot(mesh)

plt.show()
##PLOTTING
# data_2d = pvs.PVDReader(FileName=file_string)
# pvs.Show(data_2d)
# prop = pvs.GetDisplayProperties(data_2d)
# prop.Representation = "Wireframe"
# print prop.GetProperty("Representation")
# pvs.Interact(view=None)
# pvs.Render()
