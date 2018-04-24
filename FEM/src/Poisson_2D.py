import dolfin
import matplotlib.pyplot as plt
# from paraview import simple as pvs
import subprocess
import mesh_writer
import numpy as np

elec_length = 15.0e-9
elec_height = 10.0e-9
elec_spacing = 10.0e-9
boundary_x_min = 0.0e-9
boundary_x_max = (4*elec_length+4*elec_spacing)
boundary_y_min = -50.0e-9
boundary_y_max = 50.0e-9
boundary_dielectric = 10.0e-9
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
        global elec_length, elec_spacing, mid_x
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+2*(elec_length+elec_spacing), mid_x+1.999*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) or \
           (dolfin.between(x[0], (mid_x-1.999*(elec_length+elec_spacing), mid_x+0.5*elec_length-2*(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
        )

class Electrode_90(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length-(elec_length+elec_spacing), mid_x+0.5*elec_length-(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
        )

class Electrode_180(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length, mid_x+0.5*elec_length)) and \
            dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
        )

class Electrode_270(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global elec_length, elec_spacing, mid_x
        return ( \
           (dolfin.between(x[0], (mid_x-0.5*elec_length+(elec_length+elec_spacing), mid_x+0.5*elec_length+(elec_length+elec_spacing))) and \
            dolfin.between(x[1], (-0.5*elec_height, 0.5*elec_height))) \
        )

# INTERNAL BOUNDARY CONDITION
#DIELECTRIC
class Air(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        global boundary_y_max
        return (dolfin.between(x[1], (boundary_dielectric, boundary_y_max)) )

# Sub domain for Periodic boundary condition
class PeriodicBoundary(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] < dolfin.DOLFIN_EPS + boundary_x_min and x[0] > -(dolfin.DOLFIN_EPS + boundary_x_min) and on_boundary)
    def map(self, x, y):
        #map the left domain to the right domain to boundary_x_max
        y[0] = x[0] - boundary_x_max 
        y[1] = x[1]

# Create periodic boundary condition
pbc = PeriodicBoundary()

# Initialize sub-domain instances
print "Initialising subdomain instances..."
left = Left()
top = Top()
right = Right()
bottom = Bottom()
electrode_0 = Electrode_0()
electrode_90 = Electrode_90()
electrode_180 = Electrode_180()
electrode_270 = Electrode_270()
air = Air()

# Use in-house package to define mesh
print "Initializing mesh with MeshWriter..."
mw = mesh_writer.MeshWriter()
mw.resolution = boundary_x_max/10.0
mw.addOuterBound([boundary_x_min,boundary_y_min], [boundary_x_max,boundary_y_max], 1)
mw.addCrack([0.01*boundary_x_max,boundary_dielectric],[0.999*boundary_x_max,boundary_dielectric],0.05)
mw.addCrack([0.01*boundary_x_max,0.8*boundary_dielectric],[0.999*boundary_x_max,0.8*boundary_dielectric],0.1)

mw.addCrackBox([mid_x-0.5*elec_length,-0.5*elec_height],\
               [mid_x+0.5*elec_length,0.5*elec_height],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length + mid_x+0.5*elec_length)/2.0, 0.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length-(elec_length+elec_spacing),-0.5*elec_height],\
               [mid_x+0.5*elec_length-(elec_length+elec_spacing),0.5*elec_height],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length-(elec_length+elec_spacing) + mid_x+0.5*elec_length-(elec_length+elec_spacing))/2.0, 0.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length+(elec_length+elec_spacing),-0.5*elec_height],\
               [mid_x+0.5*elec_length+(elec_length+elec_spacing),0.5*elec_height],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length+(elec_length+elec_spacing) + mid_x+0.5*elec_length+(elec_length+elec_spacing))/2.0, 0.0], 1)

#Use 1.99 instead of 2.00 to prevent intersection of face with boundary face.
mw.addCrackBox([mid_x-1.999*(elec_length+elec_spacing),-0.5*elec_height],\
               [mid_x+0.5*elec_length-2*(elec_length+elec_spacing),0.5*elec_height],0.1)
mw.addPointToSurface([(mid_x-1.999*(elec_length+elec_spacing) + mid_x+0.5*elec_length-2*(elec_length+elec_spacing))/2.0, 0.0], 1)

mw.addCrackBox([mid_x-0.5*elec_length+2*(elec_length+elec_spacing),-0.5*elec_height],\
               [mid_x+1.999*(elec_length+elec_spacing),0.5*elec_height],0.1)
mw.addPointToSurface([(mid_x-0.5*elec_length+2*(elec_length+elec_spacing) + mid_x+1.999*(elec_length+elec_spacing))/2.0, 0.0], 1)

with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
subprocess.call(['gmsh -2 ../data/domain.geo -string "General.ExpertMode=1;"'], shell=True) #Expert mode to suppress warnings about fine mesh
subprocess.call(['dolfin-convert ../data/domain.msh ../data/domain.xml'], shell=True)
mesh = dolfin.Mesh('../data/domain.xml')

# Initialize mesh function for interior domains
print "Marking domains..."
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
air.mark(domains, 1)

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


# Define function space and basis functions
print "Defining function space and basis..."
V = dolfin.FunctionSpace(mesh, "CG", 3, constrained_domain=pbc)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

# Define Dirichlet boundary conditions at top, bottom, front, back
print "Defining Dirichlet boundaries..."

############# This entire part is kinda ehhhhhhhhh...
chi_si = 4.05 #eV
phi_gold = 5.1 #eV
# phi_bi = phi_gold - chi_si
phi_bi = 0
amp = 1.05
#############

bcs = [dolfin.DirichletBC(V, amp*np.sin(0)+phi_bi, boundaries, 5),
       dolfin.DirichletBC(V, amp*np.sin(np.pi/2.0)+phi_bi, boundaries, 6),
       dolfin.DirichletBC(V, amp*np.sin(np.pi)+phi_bi, boundaries, 7),
       dolfin.DirichletBC(V, amp*np.sin(3.0/2.0*np.pi)+phi_bi, boundaries, 8)]
       # dolfin.DirichletBC(V, np.abs(amp/boundary_y_max), boundaries, 2),
       # dolfin.DirichletBC(V, np.abs(amp/boundary_y_min), boundaries, 4),
       # dolfin.DirichletBC(V, 0, boundaries, 2)]
       # dolfin.DirichletBC(V, 0, boundaries, 4)]
       # dolfin.DirichletBC(V, 0, boundaries, 1),
       # dolfin.DirichletBC(V, 0, boundaries, 3)]

# Define input data
print "Defining inputs..."
EPS_0 = 8.854e-12       #Absolute permittivity, [Farad/metre]
Q_E = 1.602e-19         #Elementary charge, [Coulomb/charge]
k = 8.617e-5            #Boltzmann constant in eV/K
T = 4                   #Temperature in Kelvin
N_c = 2.82e25           #N_c for silicon in 1/m^3
N_v = 1.83E25           #N_v for silicon in 1/m^3
n_i = np.sqrt(N_c*N_v)*np.power(T,1.5)*np.exp(-1.1/(2.0*k*T))
print "n_i = ", n_i
# a0 = dolfin.Constant(11.68*EPS_0) #Permittivity, Si
# a1 = dolfin.Constant(1.0*EPS_0) #Permittivity, Air
a0 = dolfin.Constant(11.68*EPS_0) #Permittivity, Si
a1 = dolfin.Constant(1.0*EPS_0) #Permittivity, Air
g_T = dolfin.Constant("0.0")
g_B = dolfin.Constant("0.0")
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
#with T = 4 kelvin, k = 8.617e-5 ev/K, carrier charge density is ~0 inside silicon.
f0 = dolfin.Constant(n_i) #Temperature dependent carrier density
f1 = dolfin.Constant(0.0) #Charge density, Air
# h_T = dolfin.Constant(np.abs(1.0/boundary_y_max))
# h_B = dolfin.Constant(np.abs(1.0/boundary_y_min))
h_T = dolfin.Constant(0.0)
h_B = dolfin.Constant(0.0)
# h_L = dolfin.Constant(0.0)
# h_R = dolfin.Constant(0.0)

# Define new measures associated with the interior domains and
# exterior boundaries
print "Defining measures..."
dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)

# 1 left, 2 top, 3 right, 4 bottom
# Define variational form
print "Defining variational form..."
F = ( a0*dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + a1*dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dx(1) \
    # + h_L*u*v*ds(1) \
    + h_T*u*v*ds(2) \
    # + h_R*u*v*ds(3) \
    + h_B*u*v*ds(4) \
    # - g_L*v*ds(1) \
    # - g_T*v*ds(2) \
    # - g_R*v*ds(3) \
    # - g_B*v*ds(4) \
    - f0*v*dx(0) - f1*v*dx(1) )

# Separate left and right hand sides of equation
print "Separating LHS and RHS..."
a, L = dolfin.lhs(F), dolfin.rhs(F)

# Solve problem
u = dolfin.Function(V)
problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'gmres'
solver.parameters['preconditioner'] = 'ilu'
prm = dolfin.parameters.krylov_solver  # short form
prm.absolute_tolerance = 1E-25
prm.relative_tolerance = 1E-6
prm.maximum_iterations = 1000
print "Solving problem..."
solver.solve()
print "Solve finished"

#Print to PVD file
print "Printing solution to file..."
file_string = "../data/Potential.pvd"
dolfin.File(file_string) << u

u_P1 = dolfin.project(u, V)
u_nodal_values = u_P1.vector()
u_array = u_nodal_values.get_local()
coor = mesh.coordinates()
boundary_vals_x = []
boundary_vals_v = []
#gather the values at the Si-Air
for i in range((coor.size/2)):
    if coor[i][1]==boundary_dielectric:
        boundary_vals_x += [coor[i][0]]
        boundary_vals_v += [u(coor[i][0],coor[i][1])]
boundary_vals_x += boundary_vals_x.pop(1)
boundary_vals_v += boundary_vals_v.pop(1)
plt.figure()
plt.plot(boundary_vals_x, boundary_vals_v)
# V_vec = dolfin.VectorFunctionSpace(mesh, "CG", 3)
# gradu = dolfin.project(dolfin.grad(u),V_vec)
# plt.figure()
# plt.ylim(boundary_y_min, boundary_y_max)
# plt.xlim(boundary_x_min, boundary_x_max)
# dolfin.plot(gradu, title="gradient of u")
# #PLOTTING in matplotlib
plt.figure()
plt.ylim(boundary_y_min, boundary_y_max)
plt.xlim(boundary_x_min, boundary_x_max)
p = dolfin.plot(u, scalarbar=True, title="u")
plt.colorbar(p)
# # dolfin.plot(mesh)
# plt.figure()
# plt.ylim(boundary_y_min, boundary_y_max)
# plt.xlim(boundary_x_min, boundary_x_max)
# dolfin.plot(mesh, title="mesh")
plt.show()
print "Plotting in matplotlib..."

print "Ending."