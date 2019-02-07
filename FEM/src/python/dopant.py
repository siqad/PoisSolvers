import numpy as np
import matplotlib.pyplot as plt

# Physical constants
q = 1.60E-19 # Elementary charge - Coulomb
eps_0 = 8.85E-12 # Absolute permittivity - Farad / metre
k = 1.380E-23 # Boltzmann constant - metre^2 kilogram / second^2 Kelvin
T = 293 # Temperature - Kelvin

# Simulation parameters
resolution = 10
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

# Construct finite difference derivative matrix
ddx = -np.diag(np.ones(resolution)) + np.diag(np.ones(resolution-1), 1)
ddx[resolution-1][resolution-2] = -1
ddx[resolution-1][resolution-1] = 1
ddx = ddx/(x[1]-x[0])

# Obtain derivative of electron density
dndx = np.dot(ddx,n)

# Find the electric field that results in a drift current
# which exactly balances diffusion current
E = -k*T/q/n*dndx

# Resulting spatial charge density
# Can be obtained by taking another derivative of E, and scaling by eps.
rho = np.dot(ddx,E)

plt.plot(x,n)
plt.xlabel("X (um)") #set label text
plt.ylabel("Electron density (m^-3)")
plt.title('Electron density')

locs, labels = plt.xticks()
labels = []
for loc in locs:
    labels += [str(int(round(loc/1E-7, 2)))]
plt.xticks(locs, labels)

plt.savefig("n.pdf")
# print(E)
# print(rho)
# plt.plot(x,rho)
# plt.show()
# print(E)
# print(dndx)

# plt.plot(x, n_ext)
# plt.show()
# print(ddx)
