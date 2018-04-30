import dolfin
import subprocess
from paraview import simple as pvs
import numpy as np

boundary_x_min = 0.0
boundary_x_max = 1.0
boundary_y_min = 0.0
boundary_y_max = 1.0
boundary_z_min = 0.0
boundary_z_max = 1.0
mid_x = (boundary_x_min+boundary_x_max)/2.0
mid_y = (boundary_y_min+boundary_y_max)/2.0
mid_z = (boundary_z_min+boundary_z_max)/2.0
print mid_x, mid_y, mid_z
elec_length = 0.2
elec_width = 0.2
elec_depth = 0.2
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
    def inside(self, x, on_boundary):
        return (dolfin.between(x[0], (mid_x-elec_length/2.0, mid_x+elec_length/2.0)) \
            and dolfin.between(x[1], (mid_y-elec_width/2.0, mid_y+elec_width/2.0)) \
            and dolfin.between(x[2], (mid_z-elec_depth/2.0, mid_z+elec_depth/2.0)) )

# Initialize sub-domain instances
left = Left() #x
top = Top() #y
right = Right() #x
bottom = Bottom() #y
front = Front() #z
back = Back() #z
air = Air()
electrode = Electrode()
print "Create subdomains..."

# Define mesh
domain = """
//PARAMETERS
mesh_resolution = 0.05;
x_max = 1.0;
x_min = 0.0;
y_max = 1.0;
y_min = 0.0;
z_max = 1.0;
z_min = 0.0;

// POINTS
Point(1) = {x_min, y_min, z_min, mesh_resolution}; //0, 0, 0
Point(2) = {x_min, y_min, z_max, mesh_resolution}; //0, 0, 1
Point(3) = {x_min, y_max, z_min, mesh_resolution}; //0, 1, 0
Point(4) = {x_min, y_max, z_max, mesh_resolution}; //0, 1, 1
Point(5) = {x_max, y_min, z_min, mesh_resolution}; //1, 0, 0
Point(6) = {x_max, y_min, z_max, mesh_resolution}; //1, 0, 1
Point(7) = {x_max, y_max, z_min, mesh_resolution}; //1, 1, 0
Point(8) = {x_max, y_max, z_max, mesh_resolution}; //1, 1, 1

//2D ITEMS
//LINES
//x = x_min surface
Line(1) = {1,2}; //0, 0, 0 -> 0, 0, 1
Line(2) = {2,4}; //0, 0, 1 -> 0, 1, 1
Line(3) = {4,3}; //0, 1, 1 -> 0, 1, 0
Line(4) = {3,1}; //0, 1, 0 -> 0, 0, 0

//x = x_max surface
Line(5) = {5,6}; //1, 0, 0 -> 1, 0, 1
Line(6) = {6,8}; //1, 0, 1 -> 1, 1, 1
Line(7) = {8,7}; //1, 1, 1 -> 1, 1, 0
Line(8) = {7,5}; //1, 1, 0 -> 1, 0, 0

//lines from x_max to x_min
Line(9) = {1,5}; //0, 0, 0 -> 1, 0, 0
Line(10) = {2,6}; //0, 0, 1 -> 1, 0, 1
Line(11) = {4,8}; //0, 1, 1 -> 1, 1, 1
Line(12) = {3,7}; //0, 1, 0 -> 1, 1, 0

//SURFACES
//x = x_min
Line Loop(13) = {1,2,3,4};
Plane Surface(14) = {13};
//x = x_max
Line Loop(15) = {5,6,7,8};
Plane Surface(16) = {15};
//y = y_min
Line Loop(17) = {9,5,-10,-1};
Plane Surface(18) = {17};
//y = y_max
Line Loop(19) = {12,-7,-11,3};
Plane Surface(20) = {19};
//z = z_min
Line Loop(21) = {9,-8,-12,4};
Plane Surface(22) = {21};
//z = z_max
Line Loop(23) = {10,6,-11,-2};
Plane Surface(24) = {23};

Surface Loop(25) = {16,20,22,14,24,18};
//Surface Loop(25) = {14,18,22,16,20,24};
Volume(26) = {25};

//Electrode points
Point(9) = {0.4, 0.01, 0.4, mesh_resolution};
Point(10) = {0.6, 0.01, 0.4, mesh_resolution}; 
Point(11) = {0.4, 0.99, 0.4, mesh_resolution};
Point(12) = {0.6, 0.99, 0.4, mesh_resolution};

Point(13) = {0.4, 0.01, 0.6, mesh_resolution};
Point(14) = {0.6, 0.01, 0.6, mesh_resolution}; 
Point(15) = {0.4, 0.99, 0.6, mesh_resolution};
Point(16) = {0.6, 0.99, 0.6, mesh_resolution};

//low z face
Line(27) = {9,10};
Line(28) = {10,12};
Line(29) = {12,11};
Line(30) = {11,9};
Line Loop(31) = {27, 28, 29, 30};
Plane Surface(32) = {31};
Surface{32} In Volume{26};

//high z face
Line(33) = {13,14};
Line(34) = {14,16};
Line(35) = {16,15};
Line(36) = {15,13};
Line Loop(37) = {33, 34, 35, 36};
Plane Surface(38) = {37};
Surface{38} In Volume{26};

//low x face
Line(39) = {9,13};
Line(40) = {15,11};
Line Loop(41) = {39, -36, 40, 30};
Plane Surface(42) = {41};
Surface{42} In Volume{26};

//high x face
Line(43) = {10,14};
Line(44) = {16,12};
Line Loop(45) = {43, 34, 44, -28};
Plane Surface(46) = {45};
Surface{46} In Volume{26};

//low y face
Line Loop(47) = {27, 43, -33, -39};
Plane Surface(48) = {47};
Surface{48} In Volume{26};

//high y face
Line Loop(49) = {-29, -44, 35, 40};
Plane Surface(50) = {49};
Surface{50} In Volume{26};

//Si-Air Boundary
Point(17) = {0.01, 0.01, 0.8, mesh_resolution};
Point(18) = {0.01, 0.99, 0.8, mesh_resolution};
Point(19) = {0.99, 0.01, 0.99, mesh_resolution};
Point(20) = {0.99, 0.99, 0.99, mesh_resolution};
Line(51) = {17, 18};
Line(52) = {18, 20};
Line(53) = {20, 19};
Line(54) = {19, 17};
Line Loop(55) = {51, 52, 53, 54};
Plane Surface(56) = {55};
Surface{56} In Volume{26};

//Volume
Physical Volume(1) = {26};

"""

with open('../data/domain.geo', 'w') as f: f.write(domain)
subprocess.call(['gmsh -3 ../data/domain.geo'], shell=True)
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
electrode.mark(boundaries, 7)
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
f = dolfin.Constant(-1.0E10*Q_E)
print "Boundary values created"

# Define function space and basis functions
V = dolfin.FunctionSpace(mesh, "CG", 3)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)
print "Function, space, basis defined"

# Define Dirichlet boundary conditions at top, bottom, front, back
# bcs = [dolfin.DirichletBC(V, 0.0, boundaries, 1),
#        dolfin.DirichletBC(V, 0.0, boundaries, 3),
#        dolfin.DirichletBC(V, 0.0, boundaries, 5),
#        dolfin.DirichletBC(V, 0.0, boundaries, 6),
#        dolfin.DirichletBC(V, 1.0, boundaries, 7)]

bcs = [dolfin.DirichletBC(V, 5.0, boundaries, 7)]
print "Dirichlet boundaries defined"

# Define new measures associated with the interior domains and
# exterior boundaries
dx = dolfin.dx(subdomain_data=domains)
ds = dolfin.ds(subdomain_data=boundaries)
print "Measures defined"

# Define variational form
F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
    + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
    + dolfin.inner(a2*dolfin.grad(u), dolfin.grad(v))*dx(2) \
    - g_L*v*ds(1) - g_R*v*ds(3) \
    - g_T*v*ds(2) - g_Bo*v*ds(4) \
    - g_F*v*ds(5) - g_Ba*v*ds(6) \
    - f*v*dx(0) - f*v*dx(1) ) #- f*v*dolfin.dx(2) )
print "Variational Form defined"

# Separate left and right hand sides of equation
a, L = dolfin.lhs(F), dolfin.rhs(F)
print "Separated LHS and RHS"

# Solve problem
print "Initializing solver parameters..."
u = dolfin.Function(V)

problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
solver = dolfin.LinearVariationalSolver(problem)
# solver.parameters['linear_solver'] = 'cg'
solver.parameters['linear_solver'] = 'gmres'
solver.parameters['preconditioner'] = 'ilu'
cg_param = solver.parameters['krylov_solver']
cg_param['absolute_tolerance'] = 1E-7
cg_param['relative_tolerance'] = 1E-4
cg_param['maximum_iterations'] = 1000
print "Solving problem..."
solver.solve()
# dolfin.parameters['krylov_solver']['maximum_iterations'] = 10
# dolfin.solve(a == L, u, bcs, solver_parameters={'linear_solver': 'cg', 'preconditioner': 'ilu'})
# solve(a == L, u, bcs, solver_parameters={"newton_solver":{"relative_tolerance":1e-6}})
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
yslice.SliceType.Origin = [0, 0.5, 0]
yslice.SliceType.Normal = [0, 1, 0]
#create the slice
zslice = pvs.Slice(Input=data_3d, SliceType="Plane")  
zslice.SliceType.Origin = [0, 0, 0.5]
zslice.SliceType.Normal = [0, 0, 1]
pvs.Show(data_3d)
prop = pvs.GetDisplayProperties(data_3d)
prop.Representation = "Wireframe"
pvs.Show(zslice)
pvs.Show(yslice)
pvs.Interact(view=None)
pvs.Render()

print "Ending."