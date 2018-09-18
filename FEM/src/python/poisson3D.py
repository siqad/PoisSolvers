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
import matplotlib.colors as clrs
import siqadconn
import subdomains as sd
from dolfin_utils.meshconvert import meshconvert

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
# res_scale = float(sim_params["sim_resolution"])

vals = helpers.adjustBoundaries(boundary_x_min,boundary_x_max,\
                                boundary_y_min,boundary_y_max,\
                                metal_offset, metal_thickness)
boundary_x_min,boundary_x_max,boundary_y_min,boundary_y_max,boundary_z_min,boundary_z_max,boundary_dielectric = vals

print("Create mesh boundaries...")
bounds = [boundary_x_min,boundary_x_max,boundary_y_min,boundary_y_max,boundary_z_min,boundary_z_max,boundary_dielectric]
ps = ps_class.PoissonSolver(bounds)
ps.setSimParams(sim_params)
ps.setMetals(metal_offset, metal_thickness)
ps.setPaths(in_path=in_path, out_path=out_path)
ps.setResolution()
ps.createOuterBounds(resolution=1.0)
ps.addDielectricSurface(resolution=1.0) #dielectric seam
#over-extend a little, to ensure that the higher resolution appears.
ps.addDielectricField(res_in=0.25, res_out=1.0)

print("Create subdomains and fields...")
ps.setSubdomains()
ps.elec_list = elec_list
ps.elec_poly_list = elec_poly_list
ps.setElectrodeSubdomains()
ps.setElectrodePolySubdomains()
ps.setBGField(delta=1E-9)

print("Initializing mesh with GMSH...")
ps.writeGeoFile()
#Expert mode to suppress warnings about fine mesh
subprocess.call(["gmsh", "-3",
                os.path.join(ps.abs_in_dir,"domain.geo"),
                '-string', '"General.ExpertMode=1;"',
                '-string', '"Mesh.CharacteristicLengthFromPoints=0;"',
                '-string', '"Mesh.CharacteristicLengthExtendFromBoundary=0;"',
                '-string', '"Geometry.AutoCoherence=1;"'])
meshconvert.convert2xml(os.path.join(ps.abs_in_dir,"domain.msh"), os.path.join(ps.abs_in_dir,"domain.xml"))
mesh = dolfin.Mesh(os.path.join(ps.abs_in_dir,'domain.xml'))

######################MARKING BOUNDARIES
print("Marking boundaries...")
ps.markDomains(mesh)
# Initialize mesh function for boundary domains
ps.markBoundaries(mesh)

######################SETTING BOUNDARY VALUES
print("Creating boundary values...")
# Define input data
EPS_0 = 8.854E-12
Q_E = 1.6E-19
EPS_SI = dolfin.Constant(11.6*EPS_0)
EPS_AIR = dolfin.Constant(1.0*EPS_0)
f = dolfin.Constant("0.0")
# cap_matrix = []

mode = str(ps.sim_params["mode"])
ps.createNetlist()
steps = ps.getSteps()
for step in range(steps):
    print("Defining function, space, basis...")
    # Define function space and basis functions
    V = dolfin.FunctionSpace(mesh, "CG", 3)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    print("Defining Dirichlet boundaries...")
    ps.bcs = []
    # chi_si = 4.05 #eV
    # phi_gold = 5.1 #eV
    # phi_bi = phi_gold - chi_si
    ps.setElectrodePotentials(step, steps, V)

    ######################DEFINE MEASURES DX AND DS
    # Define new measures associated with the interior domains and
    # exterior boundaries
    print("Defining measures...")
    dx = dolfin.dx(subdomain_data=ps.domains)
    ds = dolfin.ds(subdomain_data=ps.boundaries)

    print("Defining variational form...")
    F = ( dolfin.inner(EPS_SI*dolfin.grad(u), dolfin.grad(v))*dx(0) \
        + dolfin.inner(EPS_AIR*dolfin.grad(u), dolfin.grad(v))*dx(1) \
        - f*v*dx(0) - f*v*dx(1) )

    boundary_component = ps.getBoundaryComponent(u, v, ds)
    F += boundary_component

    print("Separating LHS and RHS...")
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

    problem = dolfin.LinearVariationalProblem(a, L, u, ps.bcs)
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

    if mode == "cap":
        x0, x1, x2 = dolfin.MeshCoordinates(mesh)
        eps = dolfin.conditional(x2 <= 0.0, EPS_SI, EPS_AIR)
        cap_list = [0.0]*len(ps.net_list)
        for j in ps.net_list:
            for i in range(len(ps.elec_list)):
                curr_net = ps.elec_list[i].net
                dS = dolfin.Measure("dS")[ps.boundaries]
                n = dolfin.FacetNormal(mesh)
                m = dolfin.avg(dolfin.dot(eps*dolfin.grad(u), n))*dS(7+i)
                # average is used since +/- sides of facet are arbitrary
                v = dolfin.assemble(m)
                print("\int grad(u) * n ds({}) = ".format(7+i), v)
                cap_list[ps.net_list.index(curr_net)] = cap_list[ps.net_list.index(curr_net)] + v
            for i in range(len(ps.elec_poly_list)):
                curr_net = ps.elec_poly_list[i].net
                dS = dolfin.Measure("dS")[ps.boundaries]
                n = dolfin.FacetNormal(mesh)
                m = dolfin.avg(dolfin.dot(eps*dolfin.grad(u), n))*dS(7+len(ps.elec_list)+i)
                # average is used since +/- sides of facet are arbitrary
                v = dolfin.assemble(m)
                print("\int grad(u) * n ds({}) = ".format(7+len(ps.elec_list)+i), v)
                print(ps.net_list.index(curr_net))
                cap_list[ps.net_list.index(curr_net)] = cap_list[ps.net_list.index(curr_net)] + v
                print(cap_list)
        ps.cap_matrix.append(cap_list)

    # PRINT TO FILE
    print("Creating 2D data slice")
    X, Y, Z, nx, ny = ps.create2DSlice(u)
    u_old = u #Set the potential to use as initial guess

    print("Saving 2D potential data to XML")
    XYZ = []
    for i in range(nx):
        for j in range(ny):
            XYZ.append([X[i,j],Y[i,j],Z[i,j]])
    sqconn.export(potential=XYZ)
    if step == 0:
        ps.saveAxesPotential(X, Y, Z, "SiAirPlot.png")
        ps.saveGrad(X,Y,Z,0)
        ps.saveGrad(X,Y,Z,1)
    ps.savePotential(X,Y,Z,step)

ps.makeGif()
ps.getCaps()
#
# print("Ending...")
# print("Net list: ", ps.net_list)
