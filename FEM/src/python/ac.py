import numpy as np
import matplotlib.pyplot as plt
import os

class PowerEstimator():
    def __init__(self):
        self.p = np.empty((0,4), float)
        self.dir = ""
        self.num = 100

    # assume x is in polar form, with a radius and an angle
    def P2R(self, x):
        radius = x[0]
        angle = x[1]
        return radius * np.exp(1j*angle)

    # assume x is in cartesian form, with a real and imaginary component
    def R2P(self, x):
        return np.abs(x), np.angle(x)

    #sets the area of the design (in angstrom^2)
    def setArea(self, xs, ys):
        self.area = (max(xs) - min(xs)) * (max(ys) - min(ys))

    #normalizes the power in self.p to the area of the design
    def normalizePower(self):
        self.p = self.p / self.area #self.p now in watts/angstrom^2
        self.p = self.p * 1E8 * 1E8 #self.p now in watts/cm^2

    def save(self):
        print("SAVING")
        x = np.logspace(0, 20, num=self.num)
        np.save(os.path.join(self.dir,"x.npy"), x)
        os.mkdir(os.path.join(self.dir,"ac_powers"))
        for i in range(4):
            # plt.loglog(x, self.p[:,i], label="R{}".format(i))
            np.save(os.path.join(self.dir,"ac_powers","p{}.npy").format(i), self.p[:,i])

    #plot the 4 resistor power draws on loglog.
    def plot(self):
        x = np.logspace(0, 20, num=self.num)
        plt.figure()

        for i in range(4):
            plt.loglog(x, self.p[:,i], label="R{}".format(i))

        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power in Watts/cm^2")
        plt.title("Power per area vs Frequency")
        ax = plt.axes()
        legend = ax.legend(loc='lower right')
        plt.savefig("{}{}".format(self.dir+"/ac_sweep", ".pdf"))

    #setup the circuit using the resistance and capacitance values provided, and plot the result.
    def run(self, res, cap_m):
        self.R = res
        self.cap_m = cap_m
        self.createCircuit()
        self.normalizePower()
        self.plot()
        self.save()

    #set the voltages and sweep impedances over frequencies.
    def createCircuit(self):
        self.setVoltages()
        x = np.logspace(0, 20,num=self.num)
        for val in x:
            self.setFrequency(val)
            self.setImpedances()

    #Create the generator voltages phasors
    def setVoltages(self):
        angles = np.linspace(0,270,4)
        voltages = [np.array([1,np.deg2rad(a)]) for a in angles]
        self.v_gens = [self.P2R(v) for v in voltages]

    def setFrequency(self, freq):
        self.omega = 2*np.pi*freq

    #Create the impedances with the current frequency setting
    def setImpedances(self):
        Z = np.array(self.cap_m)
        #Make the diagonal the self capacitance
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i])
        #turn the capacitances into impedances
        Z = Z*self.omega
        #Morph the Z matrix into the MNA matrix
        for i in range(len(Z)):
            Z[i][i] = -sum(Z[i]) + (1/self.R[i])
        #find the voltage at each of the nodes after the resistor.
        v_nodes = np.linalg.solve(Z,np.divide(self.v_gens,self.R))
        powers = []
        for i in range(len(self.R)):
            voltage = self.v_gens[i] - v_nodes[i]
            current = voltage / self.R[i]
            #create the power from VI*
            pow = voltage * np.conj(current)
            powers.append(pow)
        self.p = np.append(self.p, np.array([powers]), axis=0)
