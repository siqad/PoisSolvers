import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import dolfin as df
from charge_density import ChargeDensity

class Dopant:

    def __init__(self):
        pass

    def plot(x, y, title, x_label, y_label, file_name):
        plt.clf()
        plt.plot(x,y)
        plt.xlabel(x_label) #set label text
        plt.ylabel(y_label)
        plt.title(title)
        plt.savefig("{}{}".format(file_name, ".pdf"))

    # def setPhysConstants(self):



# Physical constants
q = 1.60E-19 # Elementary charge - Coulomb
eps_0 = 8.85E-12 # Absolute permittivity - Farad / metre
k = 1.380E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin

# Simulation parameters
T = 293 # Temperature - Kelvin
resolution = 20
eps_r = 11.9
eps = eps_r*eps_0

ni_si = 1E10 # in cm^-3
ni_si = ni_si*1E2*1E2*1E2 # conversion to metre^-3

# Spatial extent
x_min = -1E-7
x_max = 1E-7
x = np.linspace(x_min, x_max, resolution)

# Intrinsic electron density
n = np.ones(resolution)
n_int = n*ni_si

# Dopant desity and profile
dopant = 1E19 # in cm^-3
dopant = dopant*1E2*1E2*1E2

# Scaling and offset for sigmoid
x_offset = -50E-9
x_scaling = 100
x_arg = x_scaling*(x-x_offset)/x_max

#Creating the doping profile
profile = 1 - 1/(1+np.exp(-x_arg))
# profile_pw = np.piecewise(x, [x < 5E-7, x>=5E-7], [lambda x: 1/(1+np.exp(-x_scaling*(x-x_offset)/x_max)), 0])
# profile_pw = np.piecewise(x, [x < 5E-7, x>=5E-7], [lambda x: 1/(1+np.exp(-x_scaling*(x-x_offset)/x_max)), 0])
n_ext = profile*dopant

# Total electron density (as a fuction of x)
# Assume full ionization
n = n_ext + n_int

# Obtain derivative of electron density
dndx = np.gradient(n,x)

# Find the electric field that results in a drift current
# which exactly balances diffusion current
E = -k*T/q/n*dndx

# print(n)
# Resulting spatial charge density
# Can be obtained by taking another derivative of E, and scaling by eps.
rho = np.gradient(E, x)*eps

# Plot what we have so far
# fig, ax = plt.subplots()
# plot(x, n, "Electron Density", "X (m)", "Electron density (m^-3)", "n")
# plot(x, E, "Electric Field", "X (m)", "Electric field (Vm^-1)", "E")
# plot(x, rho, "Charge Density", "X (m)", "Charge density (m^-3)", "rho")

# Create a general mesh
mesh = df.RectangleMesh(df.Point(min(x), min(x)), df.Point(max(x), max(x)), 64, 64)
# Possible functions that can live on this mesh
V = df.FunctionSpace(mesh, 'CG', 1)
#Get the parameter
F = df.interpolate(df.Expression('x[0]',  degree=1), V)
rho_exp = df.interpolate(ChargeDensity(F, data=rho,x=x, degree=1), V)

df.plot(rho_exp)
plt.show()
