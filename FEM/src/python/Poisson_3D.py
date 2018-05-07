import sys
sys.path.insert(0, sys.argv[3])
import dolfin
import subprocess
import numpy as np
import mesh_writer_3D as mw
import electrode_parser
import time
import os

elec_list, layer_props, sim_params = electrode_parser.xml_parse(sys.argv[1])
metal_thickness, metal_offset = electrode_parser.getZparams(layer_props)
[boundary_x_min, boundary_x_max], [boundary_y_min, boundary_y_max] = electrode_parser.getBB(elec_list)
res_scale = electrode_parser.getResolutionScale(sim_params)

#prevent the minimum values from being exactly 0. Problems arise when defining dielectric surface.
if boundary_x_min == 0:
    boundary_x_min -= 0.01*boundary_x_max
if boundary_y_min == 0:
    boundary_y_min -= 0.01*boundary_y_max
boundary_z_min = metal_offset*1.25
boundary_z_max = np.abs(0.5*boundary_z_min)
boundary_dielectric = 0.0 #at the surface.

# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(dolfin.SubDomain): #x_min
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], boundary_x_min)
class Right(dolfin.SubDomain): #x_max
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], boundary_x_max)
class Bottom(dolfin.SubDomain): #y_min
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], boundary_y_min)
class Top(dolfin.SubDomain): #y_max
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], boundary_y_max)
class Back(dolfin.SubDomain): #z_min
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], boundary_z_min)
class Front(dolfin.SubDomain): #z_max
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], boundary_z_max)

# INTERNAL BOUNDARY CONDITION FOR DIELECTRIC
class Air(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[2], (boundary_dielectric, boundary_z_max)))
# INTERNAL BOUNDARY CONDITION ELECTRODE
class Electrode(dolfin.SubDomain):
    def __init__(self, xs, ys, zs):
        self.xs = xs
        self.ys = ys
        self.zs = zs
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return (dolfin.between(x[0], (self.xs[0], self.xs[1])) \
            and dolfin.between(x[1], (self.ys[0], self.ys[1])) \
            and dolfin.between(x[2], (self.zs[0], self.zs[1])) )

print("Create mesh boundaries...")
mw = mw.MeshWriter()
mw.resolution = min((boundary_x_max-boundary_x_min), (boundary_y_max-boundary_y_min), (boundary_z_max-boundary_z_min))/res_scale
mw.addBox([boundary_x_min,boundary_y_min,boundary_z_min], [boundary_x_max,boundary_y_max,boundary_z_max], 1, "bound")

#dielectric seam
mw.addSurface([boundary_x_min+0.01*np.abs(boundary_x_min),boundary_y_min+0.01*np.abs(boundary_y_min),boundary_dielectric],\
              [boundary_x_max-0.01*np.abs(boundary_x_max),boundary_y_min+0.01*np.abs(boundary_y_min),boundary_dielectric],\
              [boundary_x_max-0.01*np.abs(boundary_x_max),boundary_y_max-0.01*np.abs(boundary_y_max),boundary_dielectric],\
              [boundary_x_min+0.01*np.abs(boundary_x_min),boundary_y_max-0.01*np.abs(boundary_y_max),boundary_dielectric],1.0, "seam")

#over-extend a little, to ensure that the higher resolution appears.
fields = []
fields += [mw.addBoxField(0.25, 1.0, \
          [boundary_x_min-0.01*np.abs(boundary_x_min), boundary_x_max+0.01*np.abs(boundary_x_max)], \
          [boundary_y_min-0.01*np.abs(boundary_y_min), boundary_y_max+0.01*np.abs(boundary_y_max)], \
          [boundary_dielectric-0.01*np.abs(boundary_z_max), boundary_dielectric+0.01*np.abs(boundary_z_max)])]
fields = [mw.addMinField(fields)]

# Initialize sub-domain instances
print("Create subdomains and fields...")
left = Left() #x
top = Top() #y
right = Right() #x
bottom = Bottom() #y
front = Front() #z
back = Back() #z
air = Air()
electrode = []
for i in range(len(elec_list)):
    electrode.append(Electrode([elec_list[i]['x1'], elec_list[i]['x2']], \
                        [elec_list[i]['y1'], elec_list[i]['y2']], \
                        [metal_offset, metal_offset+metal_thickness] ) )
    mw.addBox([elec_list[i]['x1'],elec_list[i]['y1'],metal_offset], \
              [elec_list[i]['x2'],elec_list[i]['y2'],metal_offset+metal_thickness], 1, "seam")

    #make resolution inside electrodes coarse
    fields += [mw.addBoxField(1.0, 0.0, \
              [elec_list[i]['x1'], elec_list[i]['x2']], \
              [elec_list[i]['y1'], elec_list[i]['y2']], \
              [metal_offset, metal_offset+metal_thickness])]
    fields = [mw.addMaxField(fields)]
    fields += [mw.addBoxField(0.25, 1.0, \
              [1.1*elec_list[i]['x1'], 1.1*elec_list[i]['x2']], \
              [1.1*elec_list[i]['y1'], 1.1*elec_list[i]['y2']], \
              [1.1*metal_offset, 1.1*(metal_offset+metal_thickness)])]
    fields = [mw.addMinField(fields)]

bg_field_ind = mw.addMeanField(fields, 1E-9)
mw.setBGField(bg_field_ind)
              
print("Initializing mesh with GMSH...")
with open(os.path.abspath(os.path.dirname(sys.argv[1])) + '/domain.geo', 'w') as f: f.write(mw.file_string)
subprocess.call(['gmsh -3 '+ os.path.abspath(os.path.dirname(sys.argv[1])) + '/domain.geo -string "General.ExpertMode=1;"'+\
                 ' -string "Mesh.CharacteristicLengthFromPoints=0;"'+\
                 ' -string "Mesh.CharacteristicLengthExtendFromBoundary=0;"'], shell=True) #Expert mode to suppress warnings about fine mesh
subprocess.call(['dolfin-convert ' + os.path.abspath(os.path.dirname(sys.argv[1])) + '/domain.msh ' + \
                os.path.abspath(os.path.dirname(sys.argv[1])) + '/domain.xml'], shell=True)
mesh = dolfin.Mesh(os.path.abspath(os.path.dirname(sys.argv[1])) + '/domain.xml')


print("Marking boundaries...")
# Initialize mesh function for interior domains
domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
domains.set_all(0)
air.mark(domains, 1)
# Initialize mesh function for boundary domains
boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
front.mark(boundaries, 5)
back.mark(boundaries, 6)
for i in range(len(elec_list)):        
    electrode[i].mark(boundaries, 7+i)

print("Creating boundary values...")
# Define input data
EPS_0 = 8.854E-12
Q_E = 1.6E-19
a0 = dolfin.Constant(11.6*EPS_0)
a1 = dolfin.Constant(1.0*EPS_0)
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
g_T = dolfin.Constant("0.0")
g_Bo = dolfin.Constant("0.0")
g_F = dolfin.Constant("0.0")
g_Ba = dolfin.Constant("0.0")
h_L = dolfin.Constant("0.0")
h_R = dolfin.Constant("0.0")
h_T = dolfin.Constant("0.0")
h_Bo = dolfin.Constant("0.0")
h_F = dolfin.Constant("0.0")
h_Ba = dolfin.Constant("0.0")
f = dolfin.Constant("0.0")

print("Defining function, space, basis...")
# Define function space and basis functions
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

print("Defining Dirichlet boundaries...")
bcs = []
chi_si = 4.05 #eV
phi_gold = 5.1 #eV
phi_bi = phi_gold - chi_si
for i in range(len(elec_list)):
    bcs.append(dolfin.DirichletBC(V, float(elec_list[i]['potential'])+ phi_bi, boundaries, 7+i))

# Define new measures associated with the interior domains and
# exterior boundaries
print("Defining measures...")
dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)

print("Defining variational form...")
# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    # - g_L*v*ds(1) - g_R*v*ds(3) \
    # - g_T*v*ds(2) - g_Bo*v*ds(4) \
    # - g_F*v*ds(5) - g_Ba*v*ds(6) \
    + h_L*u*v*ds(1) + h_R*u*v*ds(3) \
    + h_T*u*v*ds(2) + h_Bo*u*v*ds(4) \
    + h_F*u*v*ds(5) + h_Ba*u*v*ds(6) \
    - f*v*dx(0) - f*v*dx(1) )

print("Separateing LHS and RHS...")
# Separate left and right hand sides of equation
a, L = dolfin.lhs(F), dolfin.rhs(F)

# Solve problem
print("Initializing solver parameters...")
u = dolfin.Function(V)

problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'gmres'
solver.parameters['preconditioner'] = 'sor'
spec_param = solver.parameters['krylov_solver']
spec_param['absolute_tolerance'] = float(electrode_parser.getParam('max_abs_error', sim_params))
spec_param['relative_tolerance'] = float(electrode_parser.getParam('max_rel_error', sim_params))
spec_param['maximum_iterations'] = 2500

print("Solving problem...")
start = time.time()
solver.solve()
end = time.time()
print(("Solve finished in " + str(end-start) + " seconds."))

# PRINT TO FILE
print("Saving solution to .pvd file...")
file_string = os.path.abspath(os.path.dirname(sys.argv[1])) + "/Potential.pvd"
dolfin.parameters['allow_extrapolation'] = True
dolfin.File(file_string) << u

mid_x = (boundary_x_min+boundary_x_max)/2.0
mid_y = (boundary_y_min+boundary_y_max)/2.0
mid_z = (boundary_z_min+boundary_z_max)/2.0

print("Creating 2D data slice")
nx, ny = (200, 200)
x = np.linspace(boundary_x_min, boundary_x_max, nx)
y = np.linspace(boundary_y_min, boundary_y_max, ny)
X, Y = np.meshgrid(x, y)
z = np.array([u(i, j, boundary_dielectric) for j in y for i in x])
Z = z.reshape(nx, ny)

m_2_pix = 1.0E10 * electrode_parser.getPixPerAngstrom(elec_list)

print("Saving 2D potential data to file")
f = open(os.path.abspath(os.path.dirname(sys.argv[1])) + "/arrays_tmp.txt", "w")
for i in range(nx):
    for j in range(ny):
        f.write("{:06.2f} {:06.2f} {:06.2e}\n".format(X[i][j]*m_2_pix, Y[i][j]*m_2_pix, Z[i][j]))
print("Ending.")