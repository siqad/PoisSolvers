 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Main function for physics engine

import sys
import dolfin
import subprocess
import numpy as np
import mesh_writer_3D as mw
import poisson_class as ps_class
import helpers
import time
import os
import matplotlib.pyplot as plt
import siqadconn
import subdomains as sd
from dolfin_utils.meshconvert import meshconvert
from PIL import Image



######################PREAMBLE

in_path = sys.argv[1]
out_path = sys.argv[2]
m_per_A = 1.0E-10 #metres per angstrom
sqconn = siqadconn.SiQADConnector("PoisSolver", in_path, out_path)
metal_offset, metal_thickness = helpers.getMetalParams(sqconn)

elec_list = helpers.getElectrodeCollections(sqconn)
elec_poly_list = helpers.getElectrodePolyCollections(sqconn)
db_list = helpers.getDBCollections(sqconn)
sim_params = sqconn.getAllParameters()
[boundary_x_min, boundary_x_max], [boundary_y_min, boundary_y_max] = helpers.getBB(elec_list, elec_poly_list)
res_scale = float(sim_params["sim_resolution"])

vals = helpers.adjustBoundaries(boundary_x_min,boundary_x_max,\
                                boundary_y_min,boundary_y_max,\
                                metal_offset, metal_thickness)
boundary_x_min,boundary_x_max,boundary_y_min,boundary_y_max,boundary_z_min,boundary_z_max,boundary_dielectric = vals

bounds = [boundary_x_min,boundary_x_max,boundary_y_min,boundary_y_max,boundary_z_min,boundary_z_max,boundary_dielectric]

######################MESH AND SUBDOMAIN DEFINITION
print("Create mesh boundaries...")
ps = ps_class.PoissonSolver(bounds)
ps.setResolution(res_scale)
# mw = mw.MeshWriter()
# ps.mw.resolution = min((boundary_x_max-boundary_x_min), (boundary_y_max-boundary_y_min), (boundary_z_max-boundary_z_min))/res_scale

ps.createOuterBounds(resolution=1.0)
# ps.mw.addBox([boundary_x_min,boundary_y_min,boundary_z_min], [boundary_x_max,boundary_y_max,boundary_z_max], 1, "bound")

#dielectric seam
ps.addDielectricSurface(resolution=1.0)
# ps.mw.addSurface([boundary_x_min+0.01*np.abs(boundary_x_min),boundary_y_min+0.01*np.abs(boundary_y_min),boundary_dielectric],\
#               [boundary_x_max-0.01*np.abs(boundary_x_max),boundary_y_min+0.01*np.abs(boundary_y_min),boundary_dielectric],\
#               [boundary_x_max-0.01*np.abs(boundary_x_max),boundary_y_max-0.01*np.abs(boundary_y_max),boundary_dielectric],\
#               [boundary_x_min+0.01*np.abs(boundary_x_min),boundary_y_max-0.01*np.abs(boundary_y_max),boundary_dielectric],1.0, "seam")

#over-extend a little, to ensure that the higher resolution appears.
ps.addDielectricField(res_in=0.25, res_out=1.0)
# fields = []
# fields += [ps.mw.addBoxField(0.25, 1.0, \
#           [boundary_x_min-0.001*np.abs(boundary_x_min), boundary_x_max+0.001*np.abs(boundary_x_max)], \
#           [boundary_y_min-0.001*np.abs(boundary_y_min), boundary_y_max+0.001*np.abs(boundary_y_max)], \
#           [boundary_dielectric-0.05*np.abs(boundary_z_max), boundary_dielectric+0.05*np.abs(boundary_z_max)])]
# fields = [ps.mw.addMinField(fields)]

# Initialize sub-domain instances
print("Create subdomains and fields...")
left = sd.Left(boundary_x_min) #x
top = sd.Top(boundary_y_max) #y
right = sd.Right(boundary_x_max) #x
bottom = sd.Bottom(boundary_y_min) #y
front = sd.Front(boundary_z_max) #z
back = sd.Back(boundary_z_min) #z
air = sd.Air((boundary_dielectric, boundary_z_max))
electrode = []
for i in range(len(elec_list)):
    electrode.append(sd.Electrode([elec_list[i].x1, elec_list[i].x2], \
                                  [elec_list[i].y1, elec_list[i].y2], \
                                  [metal_offset, metal_offset+metal_thickness] ) )
    ps.addElectrode([elec_list[i].x1, elec_list[i].x2], \
                    [elec_list[i].y1, elec_list[i].y2], \
                    [metal_offset, metal_offset+metal_thickness], resolution=1.0)
    # ps.mw.addBox([elec_list[i].x1,elec_list[i].y1,metal_offset], \
    #           [elec_list[i].x2,elec_list[i].y2,metal_offset+metal_thickness], 1, "seam")
    #
    # #make resolution inside electrodes coarse
    # ps.fields += [ps.mw.addBoxField(1.0, 0.0, \
    #           [elec_list[i].x1, elec_list[i].x2], \
    #           [elec_list[i].y1, elec_list[i].y2], \
    #           [metal_offset, metal_offset+metal_thickness])]
    # ps.fields = [ps.mw.addMaxField(ps.fields)]
    # ps.fields += [ps.mw.addBoxField(0.1, 1.0, \
    #           [2.0*elec_list[i].x1, 2.0*elec_list[i].x2], \
    #           [2.0*elec_list[i].y1, 2.0*elec_list[i].y2], \
    #           [2.0*metal_offset, 2.0*(metal_offset+metal_thickness)])]
    # ps.fields = [ps.mw.addMinField(ps.fields)]

electrode_poly = []
for elec_poly in elec_poly_list:
    electrode_poly.append(sd.ElectrodePoly(elec_poly.vertex_list, \
        [metal_offset, metal_offset+metal_thickness]))
    ps.addElectrodePoly(elec_poly.vertex_list, [metal_offset, metal_thickness], resolution=0.1)
    # ps.mw.addPolygonVolume(elec_poly.vertex_list, [metal_offset, metal_thickness], 0.1)

ps.setBGField(delta=1E-9)
# bg_field_ind = ps.mw.addMeanField(ps.fields, 1E-9)
# ps.mw.setBGField(bg_field_ind)

######################MESHING WITH GMSH
print("Initializing mesh with GMSH...")
abs_in_dir = os.path.abspath(os.path.dirname(in_path))
with open(os.path.join(abs_in_dir, 'domain.geo'), 'w') as f: f.write(ps.mw.file_string)

#Expert mode to suppress warnings about fine mesh
subprocess.call(["gmsh", "-3",
                os.path.join(abs_in_dir,"domain.geo"),
                '-string', '"General.ExpertMode=1;"',
                '-string', '"Mesh.CharacteristicLengthFromPoints=0;"',
                '-string', '"Mesh.CharacteristicLengthExtendFromBoundary=0;"',
                '-string', '"Geometry.AutoCoherence=1;"'])

meshconvert.convert2xml(os.path.join(abs_in_dir,"domain.msh"), os.path.join(abs_in_dir,"domain.xml"))

mesh = dolfin.Mesh(os.path.join(abs_in_dir,'domain.xml'))

######################MARKING BOUNDARIES
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
for i in range(len(elec_poly_list)):
    electrode_poly[i].mark(boundaries, 7+len(elec_list)+i)


######################SETTING BOUNDARY VALUES
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
    for i in range(len(elec_poly_list)):
        elec_str = "ElectrodePoly "+str(i)+" is "
        if elec_poly_list[i].electrode_type == 1:
            elec_str += "clocked, "
            tot_phase = elec_poly_list[i].phase + step*360/steps
            potential_to_set = elec_poly_list[i].potential*np.sin( np.deg2rad(tot_phase) )
        else:
            elec_str += "fixed, "
            potential_to_set = elec_poly_list[i].potential
        if metal_offset > boundary_dielectric:
            elec_str += "and above the dielectric interface."
        else:
            potential_to_set += phi_bi
            elec_str += "and below the dielectric interface."
        print(elec_str)
        bcs.append(dolfin.DirichletBC(V, float(potential_to_set), boundaries, 7+len(elec_list)+i))
    ######################DEFINE MEASURES DX AND DS
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

    init_guess = str(sim_params["init_guess"])
    if init_guess == "prev":
        if step == 0:
            u = dolfin.Function(V)
        else:
            u = dolfin.interpolate(u_old,V)
    elif init_guess == "zero":
        u = dolfin.Function(V)

    problem = dolfin.LinearVariationalProblem(a, L, u, bcs)
    solver = dolfin.LinearVariationalSolver(problem)
    solver.parameters['linear_solver'] = 'gmres'
    solver.parameters['preconditioner'] = 'sor'
    spec_param = solver.parameters['krylov_solver']


    if init_guess == "prev":
        if step == 0:
            spec_param['nonzero_initial_guess'] = False
        else:
            spec_param['nonzero_initial_guess'] = True
    elif init_guess == "zero":
        spec_param['nonzero_initial_guess'] = False

    spec_param['absolute_tolerance'] = float(sim_params["max_abs_error"])
    spec_param['relative_tolerance'] = float(sim_params["max_rel_error"])
    spec_param['maximum_iterations'] = int(sim_params["max_linear_iters"])
    # spec_param['monitor_convergence'] = True
    print("Solving problem...")
    start = time.time()
    solver.solve()
    end = time.time()
    print(("Solve finished in " + str(end-start) + " seconds."))

    u.set_allow_extrapolation(True)
    if db_list:
        db_pots = []
        for db in db_list:
            db_pots.append([db.x, db.y, u(db.x*m_per_A,db.y*m_per_A, boundary_dielectric)])
        sqconn.export(db_pot=db_pots)


    # PRINT TO FILE
    abs_out_dir = os.path.abspath(os.path.dirname(out_path))
    depth = float(sim_params['slice_depth'])*1e-9
    print("Creating 2D data slice")
    nx = int(sim_params['image_resolution'])
    ny = nx
    x = np.linspace(boundary_x_min, boundary_x_max, nx)
    y = np.linspace(boundary_y_min, boundary_y_max, ny)
    X, Y = np.meshgrid(x, y)
    z = np.array([u(i, j, boundary_dielectric-depth) for j in y for i in x])
    Z = z.reshape(nx, ny)

    u_old = u

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
        savestring = os.path.join(abs_out_dir,'SiAirPlot.png')
        plt.savefig(savestring, bbox_inces="tight", pad_inches=0)
        plt.close(fig)
    fig = plt.figure(frameon=False)
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.pcolormesh(X,Y,Z,cmap=plt.cm.get_cmap('RdBu_r'))
    savestring = os.path.join(abs_out_dir,'SiAirBoundary{:03d}.png'.format(step))
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(savestring)
    plt.close(fig)
if mode == "clock":
    images = []
    image_files = []
    for file in os.listdir(os.path.dirname(in_path)):
        if file.startswith("SiAirBoundary"):
            image_files.append(os.path.join(os.path.dirname(in_path), file))
    image_files.sort()
    for image_name in image_files:
        images.append(Image.open(image_name))
    images[0].save(os.path.join(os.path.dirname(in_path), "SiAirBoundary.gif"),
               save_all=True,
               append_images=images[1:],
               delay=0.5,
               loop=0)
print("Ending...")
