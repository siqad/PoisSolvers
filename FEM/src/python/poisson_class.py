 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Class definition for physics engine

import numpy as np
import json
import os
import subdomains as sd
import dolfin
import helpers
import time
from dolfin_utils.meshconvert import meshconvert
from mesher import Mesher
from plotter import Plotter
from exporter import Exporter
from capacitance import CapacitanceEstimator
from resistance import ResistanceEstimator
from ac import PowerEstimator
from dopant import Dopant
import matplotlib.pyplot as plt

class PoissonSolver():
    #Constructor
    def __init__(self):

        #mesh generation
        self.mesher = Mesher()
        #plotting
        self.plotter = Plotter()

        #Resistance and capacitance tools
        self.cap = CapacitanceEstimator()
        self.res = ResistanceEstimator()
        self.ac = PowerEstimator()

        #empty lists and dictionaries
        self.net_list = []
        self.bcs = []

#Functions that users are expected to use.
    def initialize(self, json_export_path=None):
        print("Initialising PoissonSolver...")
        self.metal_params = helpers.getMetalParams(self.sqconn)
        self.elec_list = helpers.getElectrodeCollections(self.sqconn)
        self.db_list = helpers.getDBCollections(self.sqconn)
        self.initParameters(self.sqconn.getAllParameters())
        self.sim_params = self.sqconn.getAllParameters()
        self.createBoundaries()
        self.exporter = Exporter(in_path=self.sqconn.inputPath(), out_path=self.sqconn.outputPath(), json_path=json_export_path, sqconn=self.sqconn)

    def setupSim(self):
        print("Setting up simulation...")
        self.markDomains(self.mesh)
        # Initialize mesh function for boundary domains
        self.markBoundaries(self.mesh)
        self.setConstants()
        self.initDoping()
        self.setFunctionSpaces()
        # self.setChargeDensity()
        self.createNetlist()
        self.steps = self.getSteps()
        self.initRC()

    def setupDolfinSolver(self, step = None):
        print("Setting up Dolfin solver...")
        # self.setFunctionSpaces()
        self.setElectrodePotentials(step)
        self.setInitGuess(step)
        self.defineVariationalForm()
        self.setSolverParams()

    def solve(self):
        print("Solving problem...")
        start = time.time()
        self.solver.solve()
        end = time.time()
        print(("Solve finished in " + str(end-start) + " seconds."))

    def export(self, step = None):
        print("Exporting...")
        self.u_s.set_allow_extrapolation(True)
        self.exporter.exportDBs(step, self.steps, self.db_list, self.u_s, self.bounds['dielectric'])
        print("Creating 2D data slice...")
        data_2d = self.create2DSlice(self.u_s, self.slice_depth, self.img_res, self.bounds)
        data_2d_xz = self.create2DSliceXZ(self.u_s, self.img_res, self.bounds)
        print("Plotting 2D data slice...")
        self.plotter.plotPotential(step, data_2d, self.exporter.abs_out_dir)
        self.plotter.plotPotential(step, data_2d_xz, self.exporter.abs_out_dir, prefix="xz")
        print("Saving 2D potential data to XML...")
        self.exporter.exportPotentialXML(data_2d)
        #last step, finish off by creating gif and exporting DB potentials
        if step == self.steps-1 and self.mode == "clock":
            self.plotter.makeGif(self.mode, self.exporter.abs_in_dir, "SiAirBoundary")
            self.exporter.exportDBHistory()
            self.exporter.exportDBJSON()

    def loopSolve(self):
        for step in range(self.steps):
            self.setupDolfinSolver(step)
            self.solve()
            self.export(step)
            self.u_old = self.u_s # Set the potential to use as initial guess
            self.cap.getCaps(step, self.steps, self.u_s)
        self.res.getResistances(self.temp)
        self.ac.run(self.res.resistances, self.cap.cap_matrix)


#Functions that the user shouldn't have to call.
    def initParameters(self, params):
        self.bc_type = str(params["bcs"])
        self.eqn = str(params["eqn"])
        self.depletion_depth = float(params["depletion_depth"])
        self.doping_conc = float(params["doping_conc"])
        self.eps_r = float(params["eps_r_dielectric"])
        self.ground_depth = float(params["ground_depth"])
        self.img_res = int(params["image_resolution"])
        self.init_guess = str(params["init_guess"])
        self.material = str(params["material"])
        self.max_abs_error = float(params["max_abs_error"])
        self.max_rel_error = float(params["max_rel_error"])
        self.max_linear_iters = int(params["max_linear_iters"])
        self.method = str(params["method"])
        self.mode = str(params["mode"])
        self.pad = [float(params["padding_x"]), float(params["padding_y"])]
        self.preconditioner = str(params["preconditioner"])
        self.sim_res = float(params["sim_resolution"])
        self.slice_depth = float(params["slice_depth"])
        self.sim_steps = int(params["steps"])
        self.temp = float(params["temp"])

    def initRC(self):
        print("Setting RC params..")
        self.cap.mode = self.mode
        self.cap.mesh = self.mesh
        self.cap.eps_si = self.eps_si
        self.cap.eps_di = self.eps_di
        self.cap.net_list = self.net_list
        self.cap.elec_list = self.elec_list
        self.cap.boundaries = self.boundaries
        self.cap.dir = self.exporter.abs_in_dir
        self.res.mode = self.mode
        self.res.elec_list = self.elec_list
        self.res.dir = self.exporter.abs_in_dir
        self.res.setMaterialData(self.material)
        self.res.getInterpolant() #create interpolant for resistivity.
        self.ac.mode = self.mode
        self.ac.dir = self.exporter.abs_in_dir
        self.ac.setArea(self.xs_unpadded, self.ys_unpadded)

    # def setChargeDensity(self):
    #     print("Setting charge density...")
    #     if self.eqn == "laplace":
    #         self.f = dolfin.Constant("0.0")
    #     elif self.eqn == "poisson":
    #         #use equilibrium charge density
    #         self.dp.dopingCalc()
    #         rho = self.dp.getAsFunction(self.dp.rho, self.mesh, "x[2]")
    #         self.f = rho
            # dolfin.plot(rho)
            # plt.savefig(self.exporter.abs_in_dir+"/asdf.pdf")
        # elif self.eqn == "poisboltz":
        #     #use the poisson boltzmann definition for charge density.
        #     self.f = q*ni*exp(q/k/T*(u))*v*dx - q*ni*exp(-q/k/T*u)*v*dx - q*n_ext_exp*v*dx
    def initDoping(self):
        self.dp = Dopant(x_min=self.bounds['zmin'],
                    x_max=0,
                    nel=250,
                    depth=-self.depletion_depth,
                    conc=self.doping_conc,
                    temp=self.temp)

    #Produces a 2D data slice, used for getting data in correct format for plotting
    def create2DSlice(self, u, depth, resolution, bounds):
        x = np.linspace(bounds['xmin'], bounds['xmax'], resolution)
        y = np.linspace(bounds['ymin'], bounds['ymax'], resolution)
        X, Y = np.meshgrid(x, y)
        z = np.array([u(i, j, bounds['dielectric']+depth) for j in y for i in x])
        Z = z.reshape(resolution, resolution)
        return X, Y, Z, resolution, resolution

    def create2DSliceXZ(self, u, resolution, bounds):
        x = np.linspace(bounds['xmin'], bounds['xmax'], resolution)
        z = np.linspace(bounds['zmin'], bounds['zmax'], resolution)
        y = (bounds['ymin'] + bounds['ymax'])/2
        X, Z = np.meshgrid(x, z)
        v = np.array([u(i, y, k) for k in z for i in x])
        V = v.reshape(resolution, resolution)
        return X, Z, V, resolution, resolution


    def setSolverParams(self, step = None):
        if self.eqn == "poisboltz":
            print("PoisBoltz")
        else:
            print("Separating LHS and RHS...")
            # Separate left and right hand sides of equation
            self.a, self.L = dolfin.lhs(self.F), dolfin.rhs(self.F)

            self.problem = dolfin.LinearVariationalProblem(self.a, self.L, self.u_s, self.bcs)
            self.solver = dolfin.LinearVariationalSolver(self.problem)
            dolfin.parameters['form_compiler']['optimize'] = True
            self.solver.parameters['linear_solver'] = self.method
            self.solver.parameters['preconditioner'] = self.preconditioner
            spec_param = self.solver.parameters['krylov_solver']
            #These only accessible after spec_param available.
            if self.init_guess == "prev":
                if step == 0:
                    spec_param['nonzero_initial_guess'] = False
                else:
                    spec_param['nonzero_initial_guess'] = True
            elif self.init_guess == "zero":
                spec_param['nonzero_initial_guess'] = False
            spec_param['absolute_tolerance'] = self.max_abs_error
            spec_param['relative_tolerance'] = self.max_rel_error
            spec_param['maximum_iterations'] = self.max_linear_iters

    def setInitGuess(self, step):
        print("Setting initial guess...")
        if self.init_guess == "prev":
            if step == 0:
                self.u_s = dolfin.Function(self.V)
            else:
                self.u_s = dolfin.interpolate(self.u_old,self.V)
        elif self.init_guess == "zero":
            self.u_s = dolfin.Function(self.V)

    def setGroundPlane(self):
        self.bcs.append(dolfin.DirichletBC(self.V, float(0), self.boundaries, 6))

    def markDomains(self, mesh):
        # Initialize mesh function for interior domains
        print("Marking boundaries...")
        self.domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
        self.domains.set_all(0)
        # self.subdomains['air'].mark(self.domains, 1)

    def setChargeDensity(self):
        print("Setting charge density...")
        if self.eqn == "laplace":
            return dolfin.Constant("0.0")*self.v*self.dx
        elif self.eqn == "poisson":
            #use equilibrium charge density
            self.dp.dopingCalc()
            return self.dp.getAsFunction(self.dp.rho, self.mesh, "x[2]")*self.v*self.dx
            # self.f = rho
        elif self.eqn == "poisboltz":
            #use the poisson boltzmann definition for charge density.        self.ni_si = 1E10 # in cm^-3
            ni = 1E-14 #in ang^-3
            q = self.q
            k = 1.38064852E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin
            T = self.temp

            self.dp.dopingCalc()
            nd = dolfin.Expression('x[2] < depth + DOLFIN_EPS ? p1 : p2', \
               depth=self.dp.depth, p1=dolfin.Constant(self.dp.d_conc), p2=dolfin.Constant(0), degree=1, domain=self.mesh)

            # return q*ni*dolfin.exp(q/k/T*(self.u))*self.v*self.dx \
            #          - q*ni*dolfin.exp(-q/k/T*self.u)*self.v*self.dx \
            #          - q*nd*self.v*self.dx
            return q*ni*dolfin.exp(q/k/T*(self.u))*self.v*self.dx - q*nd*self.v*self.dx

    def defineVariationalForm(self):
        self.setGroundPlane()
        self.setMeasures()
        print("Defining variational form...")
        self.F = dolfin.inner(self.eps_exp*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx
        rho = self.setChargeDensity()



        # self.F -= self.f*self.v*self.dx
        self.F -= rho
            # + dolfin.inner(self.eps_di*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx(1) \
            # - self.f*self.v*self.dx
        # self.F = ( dolfin.inner(self.EPS_SI*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx(0) \
        #     + dolfin.inner(self.EPS_DIELECTRIC*dolfin.grad(self.u), dolfin.grad(self.v))*self.dx(1) \
        #     - self.f*self.v*self.dx(0) - self.f*self.v*self.dx(1) )
        self.F += self.getBoundaryComponent(self.u, self.v, self.ds)

        print("Separating LHS and RHS...")
        # Separate left and right hand sides of equation
        # self.a, self.L = dolfin.lhs(self.F), dolfin.rhs(self.F)

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
        self.eps0 = 8.854E-22 #in angstroms
        self.q = 1.6E-19
        self.eps_si = dolfin.Constant(11.6*self.eps0)
        self.eps_di = dolfin.Constant(self.eps_r*self.eps0)

        eps_si = dolfin.Constant(11.6*self.eps0)
        eps_ox = dolfin.Constant(self.eps_r*self.eps0)
        self.eps_exp = dolfin.Expression('x[2] < 0 + DOLFIN_EPS ? p1 : p2',
               p1=eps_si, p2=eps_ox, degree=1, domain=self.mesh)


    def createBoundaries(self):
        xs, ys, self.xs_unpadded, self.ys_unpadded = helpers.getBB(self.sqconn, self.pad)
        vals = helpers.adjustBoundaries(xs, ys, self.metal_params, self.ground_depth)
        self.setBounds(list(vals))

    def setConnector(self, connector):
        self.sqconn = connector

    def setBounds(self, bounds):
        keys = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'dielectric']
        self.bounds = dict(zip(keys, bounds))
        print(self.bounds)

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

        if (electrode.electrode_type == 1 and self.mode == "clock"):
            elec_str += "clocked, "
            tot_phase = -electrode.phase + step*360/steps
            # potential_to_set = electrode.potential*np.sin( np.deg2rad(tot_phase))
            phase_pot = electrode.potential*np.sin( np.deg2rad(tot_phase))
            potential_to_set = electrode.pot_offset + phase_pot
        else:
            elec_str += "fixed, "
            tot_phase = -electrode.phase
            phase_pot = electrode.potential*np.sin( np.deg2rad(tot_phase))
            potential_to_set = electrode.pot_offset + phase_pot
            # potential_to_set = electrode.potential*np.sin( np.deg2rad(tot_phase))
        # print(electrode.pot_offset)

        if self.metal_params[electrode.layer_id][0] > self.bounds['dielectric']:
            elec_str += "and above the dielectric interface."
        else:
            potential_to_set += phi_bi
            elec_str += "and below the dielectric interface."
        print(elec_str)
        return potential_to_set

    def getBoundaryComponent(self, u, v, ds):
        if self.bc_type == "robin":
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
        elif self.bc_type == "neumann":
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
        elif self.bc_type == "periodic":
            h_F = dolfin.Constant("0.0")
            h_Ba = dolfin.Constant("0.0")
            component =  h_F*u*v*ds(5)
            # component =  h_F*u*v*ds(5) + h_Ba*u*v*ds(6)
        return component

    def getFunctionSpace(self, mesh):
        if self.bc_type == "periodic":
            print(self.bounds['xmin'], self.bounds['xmax'], self.bounds['ymin'], self.bounds['ymax'])
            self.pbc = sd.PeriodicBoundary(self.bounds['xmin'], self.bounds['xmax'], self.bounds['ymin'], self.bounds['ymax'])
            return dolfin.FunctionSpace(mesh, "CG", 3, constrained_domain=self.pbc)
        else:
            return dolfin.FunctionSpace(mesh, "CG", 3)

    def getSteps(self):
        if self.mode == "standard":
            steps = 1
        elif self.mode == "clock":
            steps = self.sim_steps
        elif self.mode == "cap" or self.mode == "ac":
            steps = len(self.net_list)
        elif self.mode == "res":
            steps = 0
        return steps

    def createNetlist(self):
        for elec in self.elec_list:
            if elec.net not in self.net_list:
                self.net_list.append(elec.net)

    def setElectrodePotentials(self, step):
        print("Defining Dirichlet boundaries on electrodes...")
        self.bcs = []
        if self.mode == "standard" or self.mode == "clock":
            for electrode in self.elec_list:
                potential_to_set = self.getElecPotential(electrode, step, self.steps, self.elec_list.index(electrode))
                self.bcs.append(dolfin.DirichletBC(self.V, float(potential_to_set), self.boundaries, 7+self.elec_list.index(electrode)))
        elif self.mode == "cap" or self.mode == "ac":
            for electrode in self.elec_list:
                if self.net_list[step] == electrode.net:
                    self.bcs.append(dolfin.DirichletBC(self.V, float(1.0), self.boundaries, 7+self.elec_list.index(electrode)))
                else:
                    self.bcs.append(dolfin.DirichletBC(self.V, float(0.0), self.boundaries, 7+self.elec_list.index(electrode)))


### MESHING
    def createMesh(self):
        print("Defining mesh geometry...")
        self.initMesher()
        #subdomains is a dictionary keyed by 'left', 'top', etc. with the
        #face subdomains as values.
        #electrodes is a list of the electrode subdomains.
        self.subdomains, self.electrodes = self.mesher.createGeometry()
        #convert mesh to xml suitable for dolfin
        print("Converting mesh from .msh to .xml...")
        meshconvert.convert2xml(os.path.join(self.exporter.abs_in_dir,"domain.msh"), os.path.join(self.exporter.abs_in_dir,"domain.xml"))
        print("Importing mesh from .xml...")
        #set as mesh in model
        self.mesh = dolfin.Mesh(os.path.join(self.exporter.abs_in_dir,'domain.xml'))

    def initMesher(self):
        self.mesher.eqn = self.eqn
        self.mesher.resolution = self.sim_res
        self.mesher.dir = self.exporter.abs_in_dir
        self.mesher.bounds = self.bounds
        self.mesher.elec_list = self.elec_list
        self.mesher.metal_params = self.metal_params
