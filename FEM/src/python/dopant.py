import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import dolfin as df
from charge_density import ChargeDensity

class Dopant:

    def __init__(self):
        # pass
        self.setPhysConstants()

    def setUnits(self, units="mks"):
        self.units = units
        self.setPhysConstants(units)

    # Physical constants
    def setPhysConstants(self, units="mks"):
        self.q = 1.60E-19 # Elementary charge - Coulomb
        #For Boltzmann constant, use metre version so that k*T/q is always ~25mV
        self.k = 1.380E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin
        if units == "mks":
            self.eps_0 = 8.85E-12 # Absolute permittivity - Farad / metre
        elif units == "atomic":
            self.eps_0 = 8.85E-22 # Absolute permittivity - Farad / angstrom

    def plot(x, y, title, x_label, y_label, file_name):
        plt.clf()
        plt.plot(x,y)
        plt.xlabel(x_label) #set label text
        plt.ylabel(y_label)
        plt.title(title)
        plt.savefig("{}{}".format(file_name, ".pdf"))

    # Simulation parameters
    def setParameters(self, T=293, resolution=20, eps_r=11.9):
        self.T = T # Temperature - Kelvin
        self.resolution = resolution
        self.eps_r = eps_r
        self.eps = self.eps_r*self.eps_0

    # Assume only in one dimension
    def setBoundaries(self, min, max):
        self.min = min
        self.max = max

    def setIntrinsicDensity(self, intrinsic):
        self.ni_si = intrinsic

    def setDopantDensity(self, dopant):
        self.dopant = dopant

    def setDopantProfile(self, n_ext):
        self.n_ext = n_ext

    def getXSpace(self):
        return np.linspace(self.min, self.max, self.resolution)

    def calculateRho(self):
        x = self.getXSpace()

        # Intrinsic electron density
        n = np.ones(self.resolution)
        n_int = n*self.ni_si

        # Total electron density (as a fuction of x)
        # Assume full ionization
        self.n = self.n_ext + n_int

        # Obtain derivative of electron density
        self.dndx = np.gradient(self.n,x)
        # Find the electric field that results in a drift current
        # which exactly balances diffusion current
        self.E = -self.k*self.T/self.q/self.n*self.dndx
        print(np.trapz(self.E, x=x))

        # Resulting spatial charge density
        # Can be obtained by taking another derivative of E, and scaling by eps.
        self.rho = np.gradient(self.E, x)*self.eps

    def getRhoAsExpression(self, mesh):
        x = self.getXSpace()
        # Possible functions that can live on this mesh
        V = df.FunctionSpace(mesh, 'CG', 1)
        #Get the parameter
        F = df.Expression('x[0]',  degree=1)
        return ChargeDensity(F, data=self.rho,x=x, degree=1)

    def getRhoAsFunction(self, mesh):
        V = df.FunctionSpace(mesh, 'CG', 1)
        rho_exp = self.getRhoAsExpression(mesh)
        return df.interpolate(rho_exp, V)

    def plot(self):
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

if __name__ == "__main__":

    dp = Dopant()
    dp.setUnits("atomic")
    dp.setParameters(resolution=10000)
    # Spatial extent
    # x_min = -1E-8
    # x_max = 1E-8
    x_min = -1000
    x_max = 1000
    dp.setBoundaries(x_min, x_max)

    ni_si = 1E10 # in cm^-3
    # ni_si = ni_si*1E2*1E2*1E2 # conversion to metre^-3
    ni_si = ni_si*1E-8*1E-8*1E-8 # conversion to metre^-3
    dp.setIntrinsicDensity(ni_si)

    # Dopant desity and profile
    dopant = 1E19 # in cm^-3
    # dopant = dopant*1E2*1E2*1E2
    dopant = dopant*1E-8*1E-8*1E-8

    # Scaling and offset for sigmoid
    # x_offset = -50E-9
    # x_offset = -2.5E-9
    x_offset = -300
    x_scaling = 200
    x_arg = x_scaling*(dp.getXSpace()-x_offset)/x_max

    #Creating the doping profile
    profile = 1 - 1/(1+np.exp(-x_arg))
    n_ext = profile*dopant
    dp.setDopantProfile(n_ext)
    dp.calculateRho()
    mesh = df.RectangleMesh(df.Point(x_min, x_min), df.Point(x_max, x_max), 64, 64)
    rho_exp = dp.getRhoAsExpression(mesh)
    rho_func = dp.getRhoAsFunction(mesh)
    df.plot(rho_func)
    plt.show()
