import dolfin
import matplotlib.pyplot as plt
from paraview import simple as pvs
import subprocess
import mesh_writer
import numpy as np

elec_length = 1.0
elec_spacing = 1.5
boundary_x_min = 0.0
boundary_x_max = 10.0
boundary_y_min = 0.0
boundary_y_max = 30.0
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
# ELECTRODE_Phase_shift
class Electrode_0(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x, mid_y
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+2*(elec_length+elec_spacing), mid_x-0.01+2*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) or \
           (dolfin.between(x[0], (mid_x+0.01-2*(elec_length+elec_spacing), mid_x+0.5*elec_length-2*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) \
        )

class Electrode_90(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x, mid_y
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length-(elec_length+elec_spacing), mid_x+0.5*elec_length-(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) \
        )

class Electrode_180(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x, mid_y
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length, mid_x+0.5*elec_length)) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) \
        )

class Electrode_270(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x, mid_y
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+(elec_length+elec_spacing), mid_x+0.5*elec_length+(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (mid_y-0.5*elec_length, mid_y+0.5*elec_length))) \
        )

# INTERNAL BOUNDARY CONDITION
#DIELECTRIC
class Dielectric(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_y_max
        return (dolfin.between(x[1], (0.8*boundary_y_max, boundary_y_max)) )

# Initialize sub-domain instances
print "Initialising subdomai instances..."
left = Left()
top = Top()
right = Right()
bottom = Bottom()
electrode_0 = Electrode_0()
electrode_90 = Electrode_90()
electrode_180 = Electrode_180()
electrode_270 = Electrode_270()
dielectric = Dielectric()

# Use in-house package to define mesh
print "Initializing mesh with MeshWriter..."
mw = mesh_writer.MeshWriter()
mw.resolution = boundary_x_max/10.0
mw.addOuterBound([boundary_x_min,boundary_y_min], [boundary_x_max,boundary_y_max], 1)
mw.addCrack([0.01*boundary_x_max,0.8*boundary_y_max],[0.99*boundary_x_max,0.8*boundary_y_max],0.1)
mw.addCrack([0.01*boundary_x_max,0.7*boundary_y_max],[0.99*boundary_x_max,0.7*boundary_y_max],1.0)

mw.addCrackBox([mid_x-0.5*elec_length,mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length,mid_y+0.5*elec_length],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length + mid_x+0.5*elec_length)/2.0,\
                      (mid_y-0.5*elec_length + mid_y+0.5*elec_length)/2.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length-(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length-(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length-(elec_length+elec_spacing) + mid_x+0.5*elec_length-(elec_length+elec_spacing))/2.0,\
                      (mid_y-0.5*elec_length + mid_y+0.5*elec_length)/2.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length+(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length+(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length+(elec_length+elec_spacing) + mid_x+0.5*elec_length+(elec_length+elec_spacing))/2.0,\
                      (mid_y-0.5*elec_length + mid_y+0.5*elec_length)/2.0], 1)

mw.addCrackBox([mid_x+0.01-2*(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x+0.5*elec_length-2*(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addPointToSurface([(mid_x+0.01-2*(elec_length+elec_spacing) + mid_x+0.5*elec_length-2*(elec_length+elec_spacing))/2.0,\
                      (mid_y-0.5*elec_length + mid_y+0.5*elec_length)/2.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length+2*(elec_length+elec_spacing),mid_y-0.5*elec_length],\
               [mid_x-0.01+2*(elec_length+elec_spacing),mid_y+0.5*elec_length],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length+2*(elec_length+elec_spacing) + mid_x-0.01+2*(elec_length+elec_spacing))/2.0,\
                      (mid_y-0.5*elec_length + mid_y+0.5*elec_length)/2.0], 1)

with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
subprocess.call(['gmsh -2 ../data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
mesh = dolfin.Mesh('domain.xml')

# Initialize mesh function for interior domains
print "Marking domains..."
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
dielectric.mark(domains, 1)

# Initialize mesh function for boundary domains
print "Marking boundaries..."
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
electrode_0.mark(boundaries, 5)
electrode_90.mark(boundaries, 6)
electrode_180.mark(boundaries, 7)
electrode_270.mark(boundaries, 8)

# Define input data
print "Defining inputs..."
a0 = dolfin.Constant(1.0)
a1 = dolfin.Constant(100.0)
g_T = dolfin.Constant("-1.0")
g_B = dolfin.Constant("-1.0")
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
f = dolfin.Constant(0.05)

# Define function space and basis functions
print "Defining function space and basis..."
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

# Define Dirichlet boundary conditions at top, bottom, front, back
print "Defining Dirichlet boundaries..."
amp = 20
bcs = [dolfin.DirichletBC(V, amp*np.sin(0), boundaries, 5),\
       dolfin.DirichletBC(V, amp*np.sin(np.pi/2.0), boundaries, 6),\
       dolfin.DirichletBC(V, amp*np.sin(np.pi), boundaries, 7),\
       dolfin.DirichletBC(V, amp*np.sin(3.0/2.0*np.pi), boundaries, 8)]

# Define new measures associated with the interior domains and
# exterior boundaries
print "Defining measures..."
dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)

# 1 left, 2 top, 3 right, 4 bottom
# Define variational form
print "Defining variational form..."
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    - g_L*v*ds(1) - g_R*v*ds(3) \
    - g_T*v*ds(2) - g_B*v*ds(4) \
    - f*v*dx(0) - f*v*dx(1) )

# Separate left and right hand sides of equation
print "Separating LHS and RHS..."
a, L = dolfin.lhs(F), dolfin.rhs(F)

# Solve problem
u = dolfin.Function(V)
problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'cg'
solver.parameters['preconditioner'] = 'ilu'
cg_param = solver.parameters['krylov_solver']
cg_param['absolute_tolerance'] = 1E-7
cg_param['relative_tolerance'] = 1E-4
cg_param['maximum_iterations'] = 1000
print "Solving problem..."
solver.solve()
print "Solve finished"

#Print to PVD file
print "Printing solution to file..."
file_string = "../data/Potential.pvd"
dolfin.File(file_string) << u

#PLOTTING in matplotlib
print "Plotting in matplotlib..."
plt.figure()
dolfin.plot(u, title="u")
dolfin.plot(mesh)
plt.show()

##PLOTTING in paraview
# print "Plotting in paraview..."
# data_2d = pvs.PVDReader(FileName=file_string)
# pvs.Show(data_2d)
# prop = pvs.GetDisplayProperties(data_2d)
# prop.Representation = "Wireframe"
# print prop.GetProperty("Representation")
# pvs.Interact(view=None)
# pvs.Render()

print "Ending."