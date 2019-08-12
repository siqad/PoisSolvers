import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import dolfin as df
# from charge_density import ChargeDensity
from genexp import GenExp
import scipy.integrate as itg

class Dopant:

    def __init__(self, x_min=-1000, x_max=0, nel=100, depth=-600, conc=1e19, temp=293):
        self.x_min = x_min
        self.x_max = x_max
        self.nel = nel
        self.depth = depth
        self.conc = conc
        self.temp = temp
        self.setPhysConstants()

    def setUnits(self, units="mks"):
        self.units = units
        self.setPhysConstants(units)

    # Physical constants
    def setPhysConstants(self, units="mks"):
        self.q = 1.60217662E-19 # Elementary charge - Coulomb
        #For Boltzmann constant, use metre version so that k*T/q is always ~25mV
        self.k = 1.38064852E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin
        if units == "mks":
            self.eps_0 = 8.85E-12 # Absolute permittivity - Farad / metre
        elif units == "atomic":
            self.eps_0 = 8.85E-22 # Absolute permittivity - Farad / angstrom

    def plot(x, y, title, x_label, y_label, file_name):
        # print("Plotting")
        plt.clf()
        plt.plot(x,y)
        plt.xlabel(x_label) #set label text
        plt.ylabel(y_label)
        plt.title(title)
        # plt.savefig("{}{}".format(file_name, ".pdf"))

    # Simulation parameters
    def setBulkParameters(self, T=293, resolution=20, eps_r=11.9, sim_params=None):
        # if sim_params == None:
        #     sys.exit(0)
        self.T = T # Temperature - Kelvin
        self.resolution = resolution

    def setEps(self, eps_bulk, eps_di):
        #Assume x < 0 is silicon, x > 0 is dielectric.
        x = self.getXSpace()
        eps_di_arr = np.heaviside(x, 1)*eps_di*self.eps_0
        eps_b_arr = (1 - np.heaviside(x, 0))*eps_bulk*self.eps_0
        eps_arr = eps_di_arr + eps_b_arr
        self.eps = eps_arr

    # Assume only in one dimension
    def setBoundaries(self, min, max):
        self.min = min
        self.max = max

    def setIntrinsicDensity(self, intrinsic):
        self.ni_si = intrinsic
        n = np.heaviside(self.getXSpace()+600, 1)
        n -= np.heaviside(self.getXSpace(), 0)
        self.n_int = n*self.ni_si

    def setDopantProfile(self, n_ext):
        self.n_ext = n_ext

    def setN(self, data=None):
        #Full ionization, still intrinsic where n_ext = 0
        if data is None:
            self.n = self.n_ext + self.n_int
        #Use provided data
        else:
            self.n = data

    def getXSpace(self):
        return np.linspace(self.min, self.max, self.resolution)

    def calculateRho(self):
        x = self.getXSpace()

        # self.n = self.n_ext + self.n_int

        # Obtain derivative of electron density
        self.dndx = np.gradient(self.n,x)

        # Find the electric field that results in a drift current
        # which exactly balances diffusion current
        self.E = -self.k*self.T/self.q/self.n*self.dndx
        self.v = -itg.cumtrapz(self.E, x=self.getXSpace(), initial=0)
        # Resulting spatial charge density, goes to Poisson's Equation
        self.rho = np.gradient(self.E*self.eps, x)

    def getAsExpression(self, data, mesh, axis):
        x = self.getXSpace()
        V = df.FunctionSpace(mesh, 'CG', 1)
        F = df.Expression(axis, degree=1)
        return GenExp(F, data=data, x=x, degree=1)

    def getAsFunction(self, data, mesh, axis):
        V = df.FunctionSpace(mesh, 'CG', 1)
        data_exp = self.getAsExpression(data, mesh, axis)
        return df.interpolate(data_exp, V)

    def plotFigures(self):
        # Plot what we have so far
        fig, ax = plt.subplots()
        x = self.getXSpace()
        if self.units == "mks":
            plot(x, n, "Electron Density", "X (m)", "Electron density (m^-3)", "n")
            plot(x, E, "Electric Field", "X (m)", "Electric field (V m^-1)", "E")
            plot(x, rho, "Charge Density", "X (m)", "Charge density (m^-3)", "rho")
        elif self.units == "atomic":
            plot(x, n, "Electron Density", "X (angstrom)", "Electron density (angstrom^-3)", "n")
            plot(x, E, "Electric Field", "X (angstrom)", "Electric field (V angstrom^-1)", "E")
            plot(x, rho, "Charge Density", "X (angstrom)", "Charge density (angstrom^-3)", "rho")
        plt.ion()
        plt.show()

    def dopingCalc(self):
        self.setUnits("atomic")
        self.setBulkParameters(T=self.temp, resolution=self.nel)
        self.setBoundaries(self.x_min, self.x_max)
        self.setEps(11.9, 3.9)
        self.ni_si = 1E10 # in cm^-3
        self.ni_si = self.ni_si*1E-8*1E-8*1E-8 # conversion to ang^-3
        self.setIntrinsicDensity(self.ni_si)
        # Dopant desity and profile
        self.d_conc = self.conc # in cm^-3
        self.d_conc = self.d_conc*1E-8*1E-8*1E-8 #conversion to ang^-3
        # Scaling and offset for sigmoid
        x_offset = -50E-9
        x_offset = -2.5E-9
        x_offset = self.depth
        x_scaling = 0.1
        x_arg = x_scaling*(self.getXSpace()-x_offset)
        #Creating the doping profile
        profile = 1 - 1/(1+np.exp(-x_arg))
        n_ext = profile*self.d_conc
        self.setDopantProfile(n_ext)
        self.setN()
        self.calculateRho()

if __name__ == "__main__":

    dp = Dopant()

    dp.dopingCalc()
    plt.plot(dp.rho)
    plt.show()


    # dp.setUnits("atomic")
    # dp.setBulkParameters(resolution=10000)
    #
    # # Spatial extent
    # x_min = -1000
    # x_max = 0
    # dp.setBoundaries(x_min, x_max)
    # dp.setEps(11.9, 3.9)
    #
    # ni_si = 1E10 # in cm^-3
    # ni_si = ni_si*1E-8*1E-8*1E-8 # conversion to ang^-3
    # dp.setIntrinsicDensity(ni_si)
    #
    # # Dopant desity and profile
    # dopant = 1E19 # in cm^-3
    # dopant = dopant*1E-8*1E-8*1E-8
    #
    # # Scaling and offset for sigmoid
    # x_offset = -50E-9
    # x_offset = -2.5E-9
    # x_offset = -600
    # x_scaling = 0.1
    # x_arg = x_scaling*(dp.getXSpace()-x_offset)
    #
    # #Creating the doping profile
    # profile = 1 - 1/(1+np.exp(-x_arg))
    # n_ext = profile*dopant
    # dp.setDopantProfile(n_ext)
    # dp.calculateRho()
    # mesh = df.BoxMesh(df.Point(x_min, x_min, x_min), df.Point(x_max, x_max, x_max), 2, 2, 200)
    # rho_exp = dp.getAsExpression(dp.rho, mesh, 'x[0]')
    # rho_func = dp.getAsFunction(dp.rho, mesh, 'x[0]')
