import numpy as np
import matplotlib.pyplot as plt

def plot(x, y, title, x_label, y_label, file_name):
    plt.clf()
    plt.plot(x,y)
    plt.xlabel(x_label) #set label text
    plt.ylabel(y_label)
    plt.title(title)

    locs, labels = plt.xticks()
    labels = []
    for loc in locs:
        labels += [str(int(round(loc/1E-7, 2)))]
    plt.xticks(locs, labels)

    # plt.savefig("n.pdf")
    plt.savefig("{}{}".format(file_name, ".pdf"))

# Physical constants
q = 1.60E-19 # Elementary charge - Coulomb
eps_0 = 8.85E-12 # Absolute permittivity - Farad / metre
k = 1.380E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin
T = 293 # Temperature - Kelvin

# Simulation parameters
resolution = 200
eps_r = 11.9
eps = eps_r*eps_0

ni_si = 1E10 # in cm^-3
ni_si = ni_si*1E2*1E2*1E2 # conversion to metre^-3
# Spatial extent
x_min = -1E-6
x_max = 1E-6
x = np.linspace(x_min, x_max, resolution)

# Intrinsic electron density
n = np.ones(resolution)
n_int = n*ni_si

# Dopant desity and profile
dopant = 1E19 # in cm^-3
dopant = dopant*1E2*1E2*1E2

# Scaling and offset for sigmoid
x_offset = 0
x_scaling = 50
x_arg = x_scaling*(x-x_offset)/x_max

#Creating the doping profile
profile = 1/(1+np.exp(-x_arg))
n_ext = profile*dopant

# Total electron density (as a fuction of x)
# Assume full ionization
n = n_ext + n_int

# Obtain derivative of electron density
dndx = np.gradient(n,x)

# Find the electric field that results in a drift current
# which exactly balances diffusion current
E = -k*T/q/n*dndx

# Resulting spatial charge density
# Can be obtained by taking another derivative of E, and scaling by eps.
rho = np.gradient(E, x)*eps

plot(x, n, "Electron Density", "X (um)", "Electron density (m^-3)", "n")
plot(x, E, "Electric Field", "X (um)", "Electric field (Vm^-1)", "E")
plot(x, rho, "Charge Density", "X (um)", "Charge density (m^-3)", "rho")
