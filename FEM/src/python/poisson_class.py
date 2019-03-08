 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Class definition for physics engine

import numpy as np
import os
import subdomains as sd
import dolfin
import helpers
import time
import mesher
import plotter
import capacitance
import resistance
from dolfin_utils.meshconvert import meshconvert

class PoissonSolver():
    #Constructor
    def __init__(self):
        #mesh generation
        self.mesh = None
        self.mesher = mesher.Mesher()

        #I/O
        self.sqconn = None
        self.in_path = None
        self.out_path_path = None
        self.plotter = plotter.Plotter()
        self.db_hist = []

        #parameters gotten from sqconn
        self.sim_params = None
        self.elec_list = None
        self.metal_params = None
        self.db_list = None
        self.net_list = []

        #Resistance and capacitance tools
        self.cap = capacitance.CapacitanceEstimator()
        self.res = resistance.ResistanceEstimator()

        #calculated values and ones used during simulation
        self.cap_matrix = []
        self.steps = None
        self.bcs = []
        self.dx = None
        self.ds = None
        self.F = None
        self.a = None
        self.L = None
        self.u = None
        self.u_old = None
        self.v = None
        self.V = None

#Functions that users are expected to use.
    def initialize(self):
        print("Initialising PoissonSolver...")
        self.metal_params = helpers.getMetalParams(self.sqconn)
        self.elec_list = helpers.getElectrodeCollections(self.sqconn)
        self.db_list = helpers.getDBCollections(self.sqconn)
        self.sim_params = self.sqconn.getAllParameters()
        self.createBoundaries()
        self.setPaths(self.sqconn.inputPath(), self.sqconn.outputPath())


    def setupSim(self):
        print("Setting up simulation...")
        self.markDomains(self.mesh)
        # Initialize mesh function for boundary domains
        self.markBoundaries(self.mesh)
        self.setConstants()
        self.createNetlist()
        self.steps = self.getSteps()

    def setupDolfinSolver(self, step = None):
        print("Setting up Dolfin solver...")
        self.setFunctionSpaces()
        self.setElectrodePotentials(step)
        self.defineVariationalForm()
        self.setInitGuess(step)
        self.setSolverParams()
        self.initRC()

    def solve(self):
        print("Solving problem...")
        start = time.time()
        self.solver.solve()
        end = time.time()
        print(("Solve finished in " + str(end-start) + " seconds."))

    def export(self, step = None):
        print("Exporting...")
        self.u.set_allow_extrapolation(True)
        self.exportDBs(step)
        self.exportPotential(step)
        #last step, finish off by creating gif and getting capacitances if applicable
        mode = str(self.sim_params["mode"])
        if mode == "cap":
            self.getCaps(step)
        if step == self.steps-1 and mode == "clock":
            self.createGif()
            self.exportDBHistory()
            print(self.db_hist)

    def loopSolve(self):
        for step in range(self.steps):
            self.setupDolfinSolver(step)
            self.solve()
            self.export(step)
        mode = str(self.sim_params["mode"])
        if mode == "res":
            self.initRC()
            self.getRes()


#Functions that the user shouldn't have to call.
    def getRes(self):
        # print("@@ GET RES @@")
        self.res.createResGraph(float(self.sim_params["temp"]))

    def getCaps(self, step):
        # mode = str(self.sim_params["mode"])
        # if mode == "cap":
        self.cap.calcCaps()
        if step == self.steps-1:
            self.cap.formCapMatrix()
            # self.rc.getDelays(self.bounds, float(self.sim_params["temp"]))

    def initRC(self):
        print("Setting RC params..")
        mode = str(self.sim_params["mode"])
        if mode == "cap":
            self.cap.mesh = self.mesh
            self.cap.EPS_SI = self.EPS_SI
            self.cap.EPS_DIELECTRIC = self.EPS_DIELECTRIC
            self.cap.net_list = self.net_list
            self.cap.elec_list = self.elec_list
            self.cap.boundaries = self.boundaries
            self.cap.u = self.u
            self.cap.dir = self.abs_in_dir
        elif mode == "res":
            self.res.elec_list = self.elec_list
            self.res.dir = self.abs_in_dir

    def exportDBs(self, step):
        if self.db_list:
            db_pots = []
            for db in self.db_list:
                db_pots.append([2*np.pi*step/self.steps, db.x, db.y, self.u(db.x, db.y, self.bounds['dielectric'])])
            self.db_hist.extend(db_pots)
            # self.sqconn.export(db_pot=db_pots)

    def exportDBHistory(self):
            self.sqconn.export(db_pot=self.db_hist)


    def setSolverParams(self, step = None):
        self.problem = dolfin.LinearVariationalProblem(self.a, self.L, self.u, self.bcs)
        self.solver = dolfin.LinearVariationalSolver(self.problem)
        self.solver.parameters['linear_solver'] = 'gmres'
        self.solver.parameters['preconditioner'] = 'sor'
        spec_param = self.solver.parameters['krylov_solver']
        init_guess = str(self.sim_params["init_guess"])
        if init_guess == "prev":
            if step == 0:
                spec_param['nonzero_initial_guess'] = False
            else:
                spec_param['nonzero_initial_guess'] = True
        elif init_guess == "zero":
            spec_param['nonzero_initial_guess'] = False
        spec_param['absolute_tolerance'] = float(self.sim_params["max_abs_error"])
        spec_param['relative_tolerance'] = float(self.sim_params["max_rel_error"])
        spec_param['maximum_iterations'] = int(self.sim_params["max_linear_iters"])

    def setInitGuess(self, step):
        print("Setting initial guess...")
        init_guess = str(self.sim_params["init_guess"])
        if init_guess == "prev":
            if step == 0:
                self.u = dolfin.Function(self.V)
            else:
                self.u = dolfin.interpolate(self.u_old,self.V)
        elif init_guess == "zero":
            self.u = dolfin.Function(self.V)

    def setGroundPlane(self):
        self.bcs.append(dolfin.DirichletBC(self.V, float(0), self.boundaries, 6))

    def defineVariationalForm(self):
        self.setGroundPlane()
        self.setMeasures()
        print("Defining variational form...")
        self.F = ( dolfin.inner(self.EPS_SI*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx(0) \
            + dolfin.inner(self.EPS_DIELECTRIC*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx(1) \
            - self.f*self.v*self.dx(0) - self.f*self.v*self.dx(1) )
        self.F += self.getBoundaryComponent(self.u, self.v, self.ds)
        print("Separating LHS and RHS...")
        # Separate left and right hand sides of equation
        self.a, self.L = dolfin.lhs(self.F), dolfin.rhs(self.F)

    def setFunctionSpaces(self):
        print("Defining function, space, basis...")
        self.V = self.getFunctionSpace(self.mesh)
        self.u = dolfin.TrialFunction(self.V)
        self.v = dolfin.TestFunction(self.V)

    def setMeasures(self):
        print("Defining measures...")
        self.dx = dolfin.dx(subdomain_data=self.domains)
        self.ds = dolfin.ds(subdomain_data=self.boundaries)

    def setConstants(self):
        print("Setting constants...")
        # Define input data
        self.EPS_0 = 8.854E-22 #in angstroms
        self.Q_E = 1.6E-19
        self.EPS_SI = dolfin.Constant(11.6*self.EPS_0)
        EPS_R = float(self.sim_params["eps_r_dielectric"])
        # mode = str(self.sim_params["mode"]))
        self.EPS_DIELECTRIC = dolfin.Constant(EPS_R*self.EPS_0)
        print("Setting charge density...")
        self.f = dolfin.Constant("0.0")


    def createBoundaries(self):
        xs, ys = helpers.getBB(self.sqconn)
        ground_plane = float(self.sim_params["ground_depth"])
        vals = helpers.adjustBoundaries(xs, ys, self.metal_params, ground_plane)
        # vals = helpers.adjustBoundaries(xs, ys, self.metal_params)
        self.setBounds(list(vals))

    def setConnector(self,connector):
        self.sqconn = connector

    def setBounds(self, bounds):
        keys = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'dielectric']
        self.bounds = dict(zip(keys, bounds))
        print(self.bounds)

    def setPaths(self, in_path="", out_path=""):
        if in_path != "":
            self.in_path = in_path
            self.abs_in_dir = os.path.abspath(os.path.dirname(self.in_path))
        if out_path != "":
            self.out_path = out_path
            self.abs_out_dir = os.path.abspath(os.path.dirname(self.out_path))

    def setSimParams(self, sim_params=None):
        if sim_params:
            self.sim_params = sim_params

    def markDomains(self, mesh):
        # Initialize mesh function for interior domains
        print("Marking boundaries...")
        self.domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
        self.domains.set_all(0)
        self.subdomains['air'].mark(self.domains, 1)

    def markBoundaries(self, mesh):
        self.boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
        self.boundaries.set_all(0)
        self.subdomains['left'].mark(self.boundaries, 1)
        self.subdomains['top'].mark(self.boundaries, 2)
        self.subdomains['right'].mark(self.boundaries, 3)
        self.subdomains['bottom'].mark(self.boundaries, 4)
        self.subdomains['front'].mark(self.boundaries, 5)
        self.subdomains['back'].mark(self.boundaries, 6)

        #mark each electrode by its position in the list, plus an offset
        for electrode in self.electrodes:
            electrode.mark(self.boundaries, 7 + self.electrodes.index(electrode))

    def getElecPotential(self, electrode, step, steps, i):
        chi_si = 4.05 #eV
        phi_gold = 5.1 #eV
        phi_bi = phi_gold - chi_si
        elec_str = "Electrode "+str(i)+" is "
        mode = str(self.sim_params["mode"])

        if (electrode.electrode_type == 1 and mode == "clock"):
            elec_str += "clocked, "
            tot_phase = electrode.phase + step*360/steps
            potential_to_set = electrode.potential*np.sin( np.deg2rad(tot_phase) )
        else:
            elec_str += "fixed, "
            potential_to_set = electrode.potential

        if self.metal_params[electrode.layer_id][0] > self.bounds['dielectric']:
            elec_str += "and above the dielectric interface."
        else:
            potential_to_set += phi_bi
            elec_str += "and below the dielectric interface."
        print(elec_str)
        return potential_to_set

    def getBoundaryComponent(self, u, v, ds):
        if self.sim_params["bcs"] == "robin":
            h_L = dolfin.Constant("0.0")
            h_R = dolfin.Constant("0.0")
            h_T = dolfin.Constant("0.0")
            h_Bo = dolfin.Constant("0.0")
            h_F = dolfin.Constant("0.0")
            h_Ba = dolfin.Constant("0.0")
            component =  h_L*u*v*ds(1) + h_R*u*v*ds(3) \
                + h_T*u*v*ds(2) + h_Bo*u*v*ds(4) \
                + h_F*u*v*ds(5)
                # + h_F*u*v*ds(5) + h_Ba*u*v*ds(6)
        elif self.sim_params["bcs"] == "neumann":
            g_L = dolfin.Constant("0.0")
            g_R = dolfin.Constant("0.0")
            g_T = dolfin.Constant("0.0")
            g_Bo = dolfin.Constant("0.0")
            g_F = dolfin.Constant("0.0")
            g_Ba = dolfin.Constant("0.0")
            component = - g_L*v*ds(1) - g_R*v*ds(3) \
                 - g_T*v*ds(2) - g_Bo*v*ds(4) \
                 - g_F*v*ds(5)
                 # - g_F*v*ds(5) - g_Ba*v*ds(6)
        elif self.sim_params["bcs"] == "periodic":
            h_F = dolfin.Constant("0.0")
            h_Ba = dolfin.Constant("0.0")
            component =  h_F*u*v*ds(5)
            # component =  h_F*u*v*ds(5) + h_Ba*u*v*ds(6)
        return component

    def getFunctionSpace(self, mesh):
        if self.sim_params["bcs"] == "periodic":
            print(self.bounds['xmin'], self.bounds['xmax'], self.bounds['ymin'], self.bounds['ymax'])
            self.pbc = sd.PeriodicBoundary(self.bounds['xmin'], self.bounds['xmax'], self.bounds['ymin'], self.bounds['ymax'])
            return dolfin.FunctionSpace(mesh, "CG", 3, constrained_domain=self.pbc)
        else:
            return dolfin.FunctionSpace(mesh, "CG", 3)


    def getSteps(self):
        mode = str(self.sim_params["mode"])
        if mode == "standard":
            steps = 1
        elif mode == "clock":
            steps = int(self.sim_params["steps"])
        elif mode == "cap":
            steps = len(self.net_list)
        elif mode == "res":
            steps = 0
        return steps

    def createNetlist(self):
        for elec in self.elec_list:
            if elec.net not in self.net_list:
                self.net_list.append(elec.net)

    def setElectrodePotentials(self, step):
        print("Defining Dirichlet boundaries on electrodes...")
        self.bcs = []
        mode = str(self.sim_params["mode"])
        if mode == "standard" or mode == "clock":
            for electrode in self.elec_list:
                potential_to_set = self.getElecPotential(electrode, step, self.steps, self.elec_list.index(electrode))
                self.bcs.append(dolfin.DirichletBC(self.V, float(potential_to_set), self.boundaries, 7+self.elec_list.index(electrode)))
        elif mode == "cap":
            for electrode in self.elec_list:
                if self.net_list[step] == electrode.net:
                    self.bcs.append(dolfin.DirichletBC(self.V, float(1.0), self.boundaries, 7+self.elec_list.index(electrode)))
                else:
                    self.bcs.append(dolfin.DirichletBC(self.V, float(0.0), self.boundaries, 7+self.elec_list.index(electrode)))
        else:
            return

# MESHING

    def createMesh(self):
        print("Defining mesh geometry...")
        self.initMesher()
        #subdomains is a dictionary keyed by 'left', 'top', etc. with the
        #face subdomains as values.
        #electrodes is a list of the electrode subdomains.
        self.subdomains, self.electrodes = self.mesher.createGeometry()
        self.setMesh()

    def initMesher(self):
        self.mesher.resolution = float(self.sim_params["sim_resolution"])
        self.mesher.dir = self.abs_in_dir
        self.mesher.bounds = self.bounds
        self.mesher.elec_list = self.elec_list
        self.mesher.metal_params = self.metal_params

    def setMesh(self):
        print("Converting mesh from .msh to .xml...")
        meshconvert.convert2xml(os.path.join(self.abs_in_dir,"domain.msh"), os.path.join(self.abs_in_dir,"domain.xml"))
        print("Importing mesh from .xml...")
        self.mesh = dolfin.Mesh(os.path.join(self.abs_in_dir,'domain.xml'))

# PLOTTING

    def createGif(self):
        mode = str(self.sim_params["mode"])
        dir = self.abs_in_dir
        dir_files = os.listdir(os.path.dirname(self.in_path))
        self.plotter.makeGif(mode, dir, dir_files)

    def exportPotential(self, step = None):
        print("Creating 2D data slice...")
        depth = float(self.sim_params['slice_depth']) #in angstroms
        res = int(self.sim_params['image_resolution'])
        X, Y, Z, nx, ny = self.plotter.create2DSlice(self.u, depth, res, self.bounds)
        print("MAX: ", np.max(np.max(Z)))
        self.u_old = self.u #Set the potential to use as initial guess

        print("Saving 2D potential data to XML...")
        XYZ = []
        for i in range(nx):
            for j in range(ny):
                XYZ.append([X[i,j],Y[i,j],Z[i,j]])
        self.sqconn.export(potential=XYZ)
        if step == 0:
            plot_file_name = os.path.join(self.abs_out_dir,"SiAirPlot.png")
            self.plotter.saveAxesPotential(X, Y, Z, plot_file_name)
            grad_file_name = os.path.join(self.abs_out_dir,'grad.pdf')
            self.plotter.saveGrad(X,Y,Z,grad_file_name)
        pot_file_name = os.path.join(self.abs_out_dir,'SiAirBoundary{:03d}.png'.format(step))
        self.plotter.savePotential(X,Y,Z,step,pot_file_name)
