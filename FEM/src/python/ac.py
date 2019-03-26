import numpy as np
import matplotlib.pyplot as plt

# assume x is in polar form, with a radius and an angle
def P2R(x):
    radius = x[0]
    angle = x[1]
    return radius * np.exp(1j*angle)

# assume x is in cartesian form, with a real and imaginary component
def R2P(x):
    return np.abs(x), np.angle(x)

class PowerEstimator():
    def __init__(self):
        self.p = np.empty((0,4), float)
        self.dir = ""
        self.num = 100

    #sets the area (in angstrom^2)
    def setArea(self, xs, ys):
        self.area = (max(xs) - min(xs)) * (max(ys) - min(ys))


    def normalizePower(self):
        self.p = self.p / self.area #self.p now in watts/angstrom^2
        self.p = self.p * 1E8 * 1E8 #self.p now in watts/cm^2

    def plot(self):
        x = np.logspace(0, 20, num=self.num)
        plt.figure()
        plt.loglog(x, self.p[:,0], label="R0")
        plt.loglog(x, self.p[:,1], label="R1")
        plt.loglog(x, self.p[:,2], label="R2")
        plt.loglog(x, self.p[:,3], label="R3")

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power in Watts/cm^2")
        plt.title("Power per area vs Frequency")
        ax = plt.axes()
        legend = ax.legend(loc='lower right')
        plt.savefig("{}{}".format(self.dir+"/ac_sweep", ".pdf"))

    def setup(self, res, cap_m):
        self.R = res
        self.cap_m = cap_m
        self.createCircuit()
        self.normalizePower()
        self.plot()

    def createCircuit(self):
        self.setVoltages()

        x = np.logspace(0, 20,num=self.num)
        for val in x:
            self.setFrequency(val)
            self.setImpedances()

    def setVoltages(self):
        va = np.array([1, 0])
        vb = np.array([1, np.deg2rad(90)])
        vc = np.array([1, np.deg2rad(180)])
        vd = np.array([1, np.deg2rad(270)])
        self.v_gens = [P2R(va), P2R(vb), P2R(vc), P2R(vd)]

    def setFrequency(self, freq):
        self.omega = 2*np.pi*freq

    def setImpedances(self):
        Z = np.array(self.cap_m)
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i])
        Z = Z*self.omega
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i]) + (1/self.R[i])
        v_nodes = np.linalg.solve(Z,np.divide(self.v_gens,self.R))
        powers = []
        for i in range(len(self.R)):
            voltage = self.v_gens[i] - v_nodes[i]
            current = voltage / self.R[i]
            pow = voltage * np.conj(current)
            powers.append(pow)
        self.p = np.append(self.p, np.array([powers]), axis=0)
