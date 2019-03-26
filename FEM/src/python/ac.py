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
        # self.p = np.array([[]])

    #sets the area (in angstrom^2)
    def setArea(self, xs, ys):
        self.area = (max(xs) - min(xs)) * (max(ys) - min(ys))


    def normalizePower(self):
        self.p = self.p / self.area #self.p now in watts/angstrom^2
        self.p = self.p * 1E8 * 1E8 #self.p now in watts/cm^2

    def plot(self):
        x = np.logspace(0, 20, num=self.num)
        # self.p = self.p / self.area #self.p now in watts/angstrom^2
        # self.p = self.p * 1E8 * 1E8 #self.p now in watts/cm^2
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
        # x = np.logspace(0, 20)
        # # self.p = self.p / self.area #self.p now in watts/angstrom^2
        # # self.p = self.p * 1E8 * 1E8 #self.p now in watts/cm^2
        # plt.figure()
        # plt.loglog(x, self.p[:,0], label="R0")
        # plt.loglog(x, self.p[:,1], label="R1")
        # plt.loglog(x, self.p[:,2], label="R2")
        # plt.loglog(x, self.p[:,3], label="R3")
        #
        # plt.xlabel("Frequency (Hz)")
        # plt.ylabel("Power in Watts/cm^2")
        # plt.title("Power per area vs Frequency")
        # ax = plt.axes()
        # legend = ax.legend(loc='lower right')
        # plt.savefig("{}{}".format(self.dir+"/ac_sweep", ".pdf"))

    def createCircuit(self):
        self.setVoltages()

        x = np.logspace(0, 20,num=self.num)
        for val in x:
            self.setFrequency(val)
            self.setImpedances()
        # self.setFrequency(1E20)
        # self.setImpedances()
        # self.setFrequency(1E10)
        # self.setImpedances()
        # self.setFrequency(1)
        # self.setImpedances()

    def setVoltages(self):
        va = np.array([1, 0])
        vb = np.array([1, np.deg2rad(90)])
        vc = np.array([1, np.deg2rad(180)])
        vd = np.array([1, np.deg2rad(270)])
        self.v_gens = [P2R(va), P2R(vb), P2R(vc), P2R(vd)]
        # print(self.v_gens)

    def setFrequency(self, freq):
        self.omega = 2*np.pi*freq

    def setImpedances(self):
        # Z = np.zeros([4,4])
        Z = np.array(self.cap_m)
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i])
        Z = Z*self.omega
        # print(Z)
        # Z = np.reciprocal(self.omega*Z)
        #Z is now the mutual capacitance matrix
        # print(Z)
        # A = Z
        # print("Before diagonal change")
        # print(A)
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i]) + (1/self.R[i])
        # print("After diagonal change")
        # print(Z)
        # print(A)
        # v_nodes = np.linalg.solve(A,np.divide(self.v_gens,self.R))
        v_nodes = np.linalg.solve(Z,np.divide(self.v_gens,self.R))
        powers = []
        for i in range(len(self.R)):
            voltage = self.v_gens[i] - v_nodes[i]
            # print("Voltages at {}rad/s: ".format(self.omega), voltage)
            current = voltage / self.R[i]
            pow = voltage * np.conj(current)
            powers.append(pow)
            # print(self.p)
            # print(pow)
        # print(self.p.shape)
        # print(np.array([powers]).shape)
        # self.p = np.concatenate(self.p, powers, axis=0)
        self.p = np.append(self.p, np.array([powers]), axis=0)


        # print(self.p)

            # print(power)
            # print(power)
        # print("Voltage at nodes: ")

        # print(v_nodes)
        # print(self.cap_m)



# R = [100]*4
# va = np.array([1, 0])
# vb = np.array([1, np.deg2rad(90)])
# vc = np.array([1, np.deg2rad(180)])
# vd = np.array([1, np.deg2rad(270)])
# v_gens = [P2R(va), P2R(vb), P2R(vc), P2R(vd)]
# print(v_gens)


# C_net0 = 3.4710153285650266e-19F
# C_net1 = 3.5368579719566134e-19F
# C_net2 = 3.353951273969931e-19F
# C_net3 = 3.2749458120962368e-19F
# 1.9533199170382535e-18 -8.108209265564725e-19 -1.8728167290854811e-19 -6.081157847167302e-19
# -8.108209265564725e-19 2.1734322659530066e-18 -8.163561133609602e-19 -1.9256942883991256e-19
# -1.8728167290854811e-19 -8.163561133609602e-19 2.1730140796244824e-18 -8.33981165957981e-19
# -6.081157847167302e-19 -1.9256942883991256e-19 -8.33981165957981e-19 1.9621609607242476e-18

def mess():
    # f = 1E15
    # omega = 2*np.pi*f

    Z = np.zeros([4,4])
    Z = -np.array([[1.9533199170382535e-18, -8.108209265564725e-19, -1.8728167290854811e-19, -6.081157847167302e-19],
    [-8.108209265564725e-19, 2.1734322659530066e-18, -8.163561133609602e-19, -1.9256942883991256e-19],
    [-1.8728167290854811e-19, -8.163561133609602e-19, 2.1730140796244824e-18, -8.33981165957981e-19],
    [-6.081157847167302e-19, -1.9256942883991256e-19, -8.33981165957981e-19, 1.9621609607242476e-18]])
    for i in range(len(Z)):
        Z[i][i] = -sum(Z[i])
    #Z is now the mutual capacitance matrix
    Z = np.reciprocal(omega*Z)
    # print(Z)
    #Z is now the impedance

    # A = np.zeros([4,4])
    # A = np.reciprocal(-Z)
    A = -Z
    for i in range(len(Z)):
        A[i][i] = sum(Z[i]) + (1/R[i])
    # print(A)

    v_nodes = np.linalg.solve(A,np.divide(v_gens,R))
    print(v_nodes)
