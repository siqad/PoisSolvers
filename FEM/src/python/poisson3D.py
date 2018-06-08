import sys
import dolfin
import subprocess
import numpy as np
import mesh_writer_3D as mw
import time
import os
import matplotlib.pyplot as plt
import siqadconn
from dolfin_utils.meshconvert import meshconvert

def getBB(elec_list):
    min_x = min([a.x1 for a in elec_list])
    max_x = max([a.x2 for a in elec_list])
    min_y = min([a.y1 for a in elec_list])
    max_y = max([a.y2 for a in elec_list])

    xs = [min_x-2.0*(max_x-min_x), max_x+2.0*(max_x-min_x)]
    ys = [min_y-2.0*(max_y-min_y), max_y+2.0*(max_y-min_y)]
    return xs, ys

in_path = sys.argv[1]
out_path = sys.argv[2]
sqconn = siqadconn.SiQADConnector("PoisSolver", in_path, out_path)
for layer in sqconn.getLayers():
    if layer.name == "Metal":
        metal_offset = float(layer.zoffset)
        metal_thickness = float(layer.zheight)
        break

elec_list = []
m_per_A = 1.0E-10 #metres per angstrom
for elec in sqconn.electrodeCollection():
    elec_curr = elec
    #units given are in angstroms, convert to metres
    elec_curr.x1 *= m_per_A
    elec_curr.x2 *= m_per_A
    elec_curr.y1 *= m_per_A
    elec_curr.y2 *= m_per_A
    elec_list.append(elec_curr)

sim_params = sqconn.getAllParameters()

[boundary_x_min, boundary_x_max], [boundary_y_min, boundary_y_max] = getBB(elec_list)
res_scale = float(sim_params["sim_resolution"])

#prevent the minimum values from being exactly 0. Problems arise when defining dielectric surface.
if boundary_x_min == 0:
    boundary_x_min -= 0.01*boundary_x_max
if boundary_y_min == 0:
    boundary_y_min -= 0.01*boundary_y_max
boundary_z_min = -np.abs(metal_offset)*10.0
boundary_z_max = -boundary_z_min
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
          [boundary_x_min-0.001*np.abs(boundary_x_min), boundary_x_max+0.001*np.abs(boundary_x_max)], \
          [boundary_y_min-0.001*np.abs(boundary_y_min), boundary_y_max+0.001*np.abs(boundary_y_max)], \
          [boundary_dielectric-0.05*np.abs(boundary_z_max), boundary_dielectric+0.05*np.abs(boundary_z_max)])]
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
    electrode.append(Electrode([elec_list[i].x1, elec_list[i].x2], \
                        [elec_list[i].y1, elec_list[i].y2], \
                        [metal_offset, metal_offset+metal_thickness] ) )
    mw.addBox([elec_list[i].x1,elec_list[i].y1,metal_offset], \
              [elec_list[i].x2,elec_list[i].y2,metal_offset+metal_thickness], 1, "seam")

    #make resolution inside electrodes coarse
    fields += [mw.addBoxField(1.0, 0.0, \
              [elec_list[i].x1, elec_list[i].x2], \
              [elec_list[i].y1, elec_list[i].y2], \
              [metal_offset, metal_offset+metal_thickness])]
    fields = [mw.addMaxField(fields)]
    fields += [mw.addBoxField(0.25, 1.0, \
              [1.1*elec_list[i].x1, 1.1*elec_list[i].x2], \
              [1.1*elec_list[i].y1, 1.1*elec_list[i].y2], \
              [1.1*metal_offset, 1.1*(metal_offset+metal_thickness)])]
    fields = [mw.addMinField(fields)]

bg_field_ind = mw.addMeanField(fields, 1E-9)
mw.setBGField(bg_field_ind)

print("Initializing mesh with GMSH...")
abs_in_dir = os.path.abspath(os.path.dirname(in_path))
with open(os.path.join(abs_in_dir, 'domain.geo'), 'w') as f: f.write(mw.file_string)

#Expert mode to suppress warnings about fine mesh
subprocess.call(["gmsh", "-3",
                os.path.join(abs_in_dir,"domain.geo"),
                '-string', '"General.ExpertMode=1;"',
                '-string', '"Mesh.CharacteristicLengthFromPoints=0;"',
                '-string', '"Mesh.CharacteristicLengthExtendFromBoundary=0;"'])

meshconvert.convert2xml(os.path.join(abs_in_dir,"domain.msh"), os.path.join(abs_in_dir,"domain.xml"))

mesh = dolfin.Mesh(os.path.join(abs_in_dir,'domain.xml'))

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
f = dolfin.Constant("0.0")

mode = str(sim_params["mode"])
if mode == "standard":
    steps = 1
elif mode == "clock":
    steps = int(sim_params["steps"])
# print(elec_list)
for step in range(steps):
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
        elec_str = "Electrode "+str(i)+" is "
        if elec_list[i].electrode_type == 1:
            elec_str += "clocked, "
            tot_phase = elec_list[i].phase + step*360/steps
            potential_to_set = elec_list[i].potential*np.sin( np.deg2rad(tot_phase) )
        else:
            elec_str += "fixed, "
            potential_to_set = elec_list[i].potential
        if metal_offset > boundary_dielectric:
            elec_str += "and above the dielectric interface."
        else:
            potential_to_set += phi_bi
            elec_str += "and below the dielectric interface."
        print(elec_str)
        bcs.append(dolfin.DirichletBC(V, float(potential_to_set), boundaries, 7+i))

    # Define new measures associated with the interior domains and
    # exterior boundaries
    print("Defining measures...")
    dx = dolfin.dx(subdomain_data=domains)
    ds = dolfin.ds(subdomain_data=boundaries)

    print("Defining variational form...")
    F = ( dolfin.inner(a0*dolfin.grad(u), dolfin.grad(v))*dx(0) \
        + dolfin.inner(a1*dolfin.grad(u), dolfin.grad(v))*dx(1) \
        - f*v*dx(0) - f*v*dx(1) )

    if sim_params["bcs"] == "robin":
        h_L = dolfin.Constant("0.0")
        h_R = dolfin.Constant("0.0")
        h_T = dolfin.Constant("0.0")
        h_Bo = dolfin.Constant("0.0")
        h_F = dolfin.Constant("0.0")
        h_Ba = dolfin.Constant("0.0")
        F +=  h_L*u*v*ds(1) + h_R*u*v*ds(3) \
            + h_T*u*v*ds(2) + h_Bo*u*v*ds(4) \
            + h_F*u*v*ds(5) + h_Ba*u*v*ds(6)
    elif sim_params["bcs"] == "neumann":
        g_L = dolfin.Constant("0.0")
        g_R = dolfin.Constant("0.0")
        g_T = dolfin.Constant("0.0")
        g_Bo = dolfin.Constant("0.0")
        g_F = dolfin.Constant("0.0")
        g_Ba = dolfin.Constant("0.0")
        F += - g_L*v*ds(1) - g_R*v*ds(3) \
             - g_T*v*ds(2) - g_Bo*v*ds(4) \
             - g_F*v*ds(5) - g_Ba*v*ds(6)

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
    spec_param['absolute_tolerance'] = float(sim_params["max_abs_error"])
    spec_param['relative_tolerance'] = float(sim_params["max_rel_error"])
    spec_param['maximum_iterations'] = int(sim_params["max_linear_iters"])

    print("Solving problem...")
    start = time.time()
    solver.solve()
    end = time.time()
    print(("Solve finished in " + str(end-start) + " seconds."))

    # PRINT TO FILE
    abs_out_dir = os.path.abspath(os.path.dirname(out_path))

    print("Creating 2D data slice")
    nx = int(sim_params['image_resolution'])
    ny = nx
    x = np.linspace(boundary_x_min, boundary_x_max, nx)
    y = np.linspace(boundary_y_min, boundary_y_max, ny)
    X, Y = np.meshgrid(x, y)
    z = np.array([u(i, j, boundary_dielectric) for j in y for i in x])
    Z = z.reshape(nx, ny)

    print("Saving 2D potential data to XML")
    XYZ = []
    for i in range(nx):
        for j in range(ny):
            XYZ.append([X[i,j],Y[i,j],Z[i,j]])
    sqconn.export(potential=XYZ)
    if step == 0:
        fig = plt.figure()
        plt.gca().invert_yaxis()
        plt.pcolormesh(X,Y,Z,cmap=plt.cm.get_cmap('RdBu_r'))
        cbar = plt.colorbar()
        cbar.set_label("Potential (V)")
        locs, labels = plt.yticks()
        labels = []
        for loc in locs:
            labels += [str(round(loc*1e9, 2))]
        plt.yticks(locs, labels)
        locs, labels = plt.xticks()
        labels = []
        for loc in locs:
            labels += [str(round(loc*1e9, 2))]
        plt.xticks(locs, labels)
        plt.xlabel("X (nm)")
        plt.ylabel("Y (nm)")
        savestring = os.path.join(abs_out_dir,'SiAirBoundary.png')
        plt.savefig(savestring, bbox_inces="tight", pad_inches=0)
        plt.close(fig)
    fig = plt.figure(frameon=False)
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.pcolormesh(X,Y,Z,cmap=plt.cm.get_cmap('RdBu_r'))
    savestring = os.path.join(abs_out_dir,'SiAirBoundary{}.png'.format(step))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(savestring)
    plt.close(fig)
print("Ending.")
