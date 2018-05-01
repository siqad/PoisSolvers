import dolfin
import subprocess
from paraview import simple as pvs
import numpy as np
import mesh_writer_3D as mw

boundary_x_min = 0.0
boundary_x_max = 2.0
boundary_y_min = 0.0
boundary_y_max = 2.0
boundary_z_min = 0.0
boundary_z_max = 2.0
mid_x = (boundary_x_min+boundary_x_max)/2.0
mid_y = (boundary_y_min+boundary_y_max)/2.0
mid_z = (boundary_z_min+boundary_z_max)/2.0
elec_length = 0.1
elec_width = 0.1
elec_depth = 0.1
offset = 0.3
boundary_dielectric = 0.8
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

# INTERNAL BOUNDARY CONDITION
# DIELECTRIC
class Air(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (dolfin.between(x[2], (boundary_dielectric, boundary_z_max)))

# INTERNAL BOUNDARY CONDITION
# ELECTRODE
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

mw = mw.MeshWriter()

mw.resolution = min((boundary_x_max-boundary_x_min)/10.0, (boundary_y_max-boundary_y_min)/10.0, (boundary_z_max-boundary_z_min)/20.0)
mw.addBox([boundary_x_min,boundary_y_min,boundary_z_min], [boundary_x_max,boundary_y_max,boundary_z_max], 1, "bound")
fields = []
mw.addSurface([0.01,0.01,0.8],[0.99,0.01,0.8],[0.99,0.99,0.8],[0.01,0.99,0.8],1, "seam")
fields += [mw.addTHField(0.5, 1, 0.01, 0.1)]
mw.addSurface([0.0,0.0,0.5],[1.0,0.0,0.5],[1.0,1.0,0.5],[0.0,1.0,0.5],1, "bound")
fields += [mw.addTHField(0.5, 1, 0.1, 0.3)]
mw.addMinField(fields)

# Initialize sub-domain instances
left = Left() #x
top = Top() #y
right = Right() #x
bottom = Bottom() #y
front = Front() #z
back = Back() #z
air = Air()
electrode = []
for i in range(4):
    electrode.append(Electrode([i*offset+0.1, i*offset+0.1+elec_length], \
                        [mid_y-elec_width/2.0, mid_y+elec_width/2.0], \
                        [mid_z-elec_depth/2.0, mid_z+elec_depth/2.0] ) )
    mw.addBox([i*offset+0.1,mid_y-elec_width/2.0,mid_z-elec_depth/2.0], \
              [i*offset+0.1+elec_length,mid_y+elec_width/2.0,mid_z+elec_depth/2.0], 1, "seam")

print "Create subdomains..."

with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
subprocess.call(['gmsh -3 ../data/domain.geo -string "General.ExpertMode=1;"'+\
                 ' -string "Mesh.CharacteristicLengthFromPoints=0;"'+\
                 ' -string "Mesh.CharacteristicLengthExtendFromBoundary=0;"'], shell=True) #Expert mode to suppress warnings about fine mesh
# subprocess.call(['gmsh -3 ../data/domain.geo -string "Mesh.CharacteristicLengthFromPoints=0;"'], shell=True)
subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
mesh = dolfin.Mesh('domain.xml')
print "Mesh initialized"

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
for i in range(4):        
    electrode[i].mark(boundaries, 7+i)
print "Boundaries marked"

# Define input data
EPS_0 = 8.854E-12
Q_E = 1.6E-19
a0 = dolfin.Constant(11.6*EPS_0)
a1 = dolfin.Constant(1.0*EPS_0)
a2 = dolfin.Constant(1000000*EPS_0)
g_L = dolfin.Constant("0.0")
g_R = dolfin.Constant("0.0")
g_T = dolfin.Constant("0.0")
g_Bo = dolfin.Constant("0.0")
g_F = dolfin.Constant("0.0")
g_Ba = dolfin.Constant("0.0")
h_T = dolfin.Constant(0.0)
h_B = dolfin.Constant(0.0)
f = dolfin.Constant(-1.0E10*Q_E)
print "Boundary values created"

# Define function space and basis functions
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)
print "Function, space, basis defined"

bcs = [dolfin.DirichletBC(V, 5.0, boundaries, 7),
       dolfin.DirichletBC(V, 5.0, boundaries, 8),
       dolfin.DirichletBC(V, 5.0, boundaries, 9),
       dolfin.DirichletBC(V, 5.0, boundaries, 10)]
print "Dirichlet boundaries defined"

# Define new measures associated with the interior domains and
# exterior boundaries
dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)
print "Measures defined"

# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    - g_L*v*ds(1) - g_R*v*ds(3) \
    # - g_T*v*ds(2) - g_Bo*v*ds(4) \
    - g_F*v*ds(5) - g_Ba*v*ds(6) \
    + h_T*u*v*ds(2) + h_B*u*v*ds(4) \
    - f*v*dx(0) - f*v*dx(1) )
print "Variational Form defined"

# Separate left and right hand sides of equation
a, L = dolfin.lhs(F), dolfin.rhs(F)
print "Separated LHS and RHS"

# Solve problem
print "Initializing solver parameters..."
u = dolfin.Function(V)

problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'cg'
solver.parameters['preconditioner'] = 'ilu'
cg_param = solver.parameters['krylov_solver']
cg_param['absolute_tolerance'] = 1E-7
cg_param['relative_tolerance'] = 1E-4
cg_param['maximum_iterations'] = 2500
print "Solving problem..."
solver.solve()
print "Solve finished"

# PRINT TO FILE
print "Saving solution to .pvd file..."
file_string = "../data/Potential.pvd"
dolfin.File(file_string) << u
print "Saving to .csv file..."
reader = pvs.OpenDataFile("../data/Potential.pvd")
writer = pvs.CreateWriter("../data/Potential.csv", reader)
writer.WriteAllTimeSteps = 1
writer.FieldAssociation = "Points"
writer.UpdatePipeline()
print "Solution saved."

##PLOTTING
print "Starting paraview..."
data_3d = pvs.PVDReader(FileName=file_string)
#create the slice
yslice = pvs.Slice(Input=data_3d, SliceType="Plane")  
yslice.SliceType.Origin = [0, mid_y, 0]
yslice.SliceType.Normal = [0, 1, 0]
#create the slice
zslice = pvs.Slice(Input=data_3d, SliceType="Plane")  
zslice.SliceType.Origin = [0, 0, mid_z]
zslice.SliceType.Normal = [0, 0, 1]
pvs.Show(data_3d)
prop = pvs.GetDisplayProperties(data_3d)
prop.Representation = "Wireframe"
pvs.Show(zslice)
pvs.Show(yslice)
pvs.Interact(view=None)
pvs.Render()

print "Ending."