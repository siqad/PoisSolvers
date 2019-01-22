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
#Get the I/O locations
in_path = sys.argv[1]
out_path = sys.argv[2]

#Instatiate an empty PoissonSolver object
ps = ps_class.PoissonSolver()

#Create a SiQADConnector object for I/O between the tool and the simulation engine
sqconn = siqadconn.SiQADConnector("PoisSolver", in_path, out_path)

#...and give it to PoissonSolver
ps.setConnector(sqconn)

#Initialize a bunch of stuff in PS now that we have the connector,
#including simulation parameters, defining mesh geometry, and I/O paths.
#Users should modify existing sim_params after initialize(), but before createMesh()
#to ensure their modications are reflected.
ps.initialize()

#Kick-start the mesh definition process. Defines the boundaries, and inserts
#electrodes surfaces into the mesh, adding coarseness guides as well.
#Finishes by writing the mesh definition to a .geo file, and converting it to
#a dolfin compatible .xml format.
ps.createMesh()

#Marks the subdomains and boundaries as top, left, bottom, right, dielectric, silicon, metal, etc.
#Also sets constants (elementary charge, permittivity, etc.)
ps.setupSim()

#At this point, users can override the spatial charge density (defaults to 0 everywhere).
#Using a dolfin expression could work, for symbolic definition of the charge density.
# ps.f = dolfin.Constant("0.0")


# mode = str(ps.sim_params["mode"])
ps.createNetlist()
steps = ps.getSteps()
for step in range(steps):
    print("Defining function, space, basis...")
    # Define function space and basis functions
    V = ps.getFunctionSpace(ps.mesh)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)

    print("Defining Dirichlet boundaries...")
    ps.bcs = []
    ps.setElectrodePotentials(step, steps, V)

    ######################DEFINE MEASURES DX AND DS
    # Define new measures associated with the interior domains and
    # exterior boundaries
    print("Defining measures...")
    dx = dolfin.dx(subdomain_data=ps.domains)
    ds = dolfin.ds(subdomain_data=ps.boundaries)

    print("Defining variational form...")
    F = ( dolfin.inner(ps.EPS_SI*dolfin.grad(u), dolfin.grad(v))*dx(0) \
        + dolfin.inner(ps.EPS_AIR*dolfin.grad(u), dolfin.grad(v))*dx(1) \
        - ps.f*v*dx(0) - ps.f*v*dx(1) )

    boundary_component = ps.getBoundaryComponent(u, v, ds)
    F += boundary_component

    print("Separating LHS and RHS...")
    # Separate left and right hand sides of equation
    a, L = dolfin.lhs(F), dolfin.rhs(F)

    # Solve problem
    print("Initializing solver parameters...")

    init_guess = str(ps.sim_params["init_guess"])
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

    spec_param['absolute_tolerance'] = float(ps.sim_params["max_abs_error"])
    spec_param['relative_tolerance'] = float(ps.sim_params["max_rel_error"])
    spec_param['maximum_iterations'] = int(ps.sim_params["max_linear_iters"])
    # spec_param['monitor_convergence'] = True
    print("Solving problem...")
    start = time.time()
    solver.solve()
    end = time.time()
    print(("Solve finished in " + str(end-start) + " seconds."))

    u.set_allow_extrapolation(True)
    if ps.db_list:
        db_pots = []
        for db in ps.db_list:
            db_pots.append([db.x, db.y, u(db.x, db.y, ps.bounds['dielectric'])])
            # db_pots.append([db.x, db.y, u(db.x, db.y, boundary_dielectric)])
        sqconn.export(db_pot=db_pots)

    ps.calcCaps(u, ps.mesh, ps.EPS_SI, ps.EPS_AIR)

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
