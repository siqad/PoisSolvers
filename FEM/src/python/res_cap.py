import numpy as np
import dolfin
class ResCap():
    def __init__(self, t = None, r = None):
        self.createData(t, r)
        self.getInterpolant()
        self.cap_matrix = []

    #Sets two lists as the dataset, combining them into a dictionary.
    #If no arguments are provided, uses cryogenic data for cobalt.
    def createData(self,t = None, r = None):
        #There was data given, use that instead
        if t != None and r != None:
            self.temps = t
            self.resistivities = r
        else:
        #Otherwise use data from White and Woods, 1959
            self.temps = [10,15,20,25,30, \
                     40,50,60,70,80, \
                     90,100,120,140,160, \
                     180,200,220,250,273,295] #Kelvin
            self.resistivities = [0.0901,0.0917,0.0956,0.103,0.116, \
                             0.161,0.234,0.339,0.469,0.629, \
                             0.809,0.999,1.409,1.869,2.349, \
                             2.839,3.319,3.809,4.589,5.239,5.889] #E-6 ohm cm
    #Get the polynomial coefficients of the interpolant.
    def getInterpolant(self, deg = None):
        if deg == None:
            deg = len(self.temps) - 1
        self.poly_coeffs = np.polyfit(np.array(self.temps),np.array(self.resistivities), deg)

    #Approximate
    def approxRes(self, x):
        return np.polyval(self.poly_coeffs, x)

    def formCapMatrix(self):
        matrix_string = ""
        cap_matrix = np.array(self.cap_matrix)
        cap_matrix = (np.transpose(cap_matrix) + cap_matrix)/2.0
        for i in range(len(cap_matrix)):
            tot_cap = 0
            for cap in cap_matrix[i]:
                tot_cap = tot_cap+cap
                matrix_string = matrix_string + "{} ".format(cap)
            print("C_net{} = {}F".format(self.net_list[i],tot_cap))
            matrix_string = matrix_string + "\n"
        print(matrix_string)

    def calcCaps(self):
        x0, x1, x2 = dolfin.MeshCoordinates(self.mesh)
        eps = dolfin.conditional(x2 <= 0.0, self.EPS_SI, self.EPS_AIR)
        cap_list = [0.0]*len(self.net_list)
        for electrode in self.elec_list:
            curr_net = electrode.net
            dS = dolfin.Measure("dS")[self.boundaries]
            n = dolfin.FacetNormal(self.mesh)
            m = dolfin.avg(dolfin.dot(eps*dolfin.grad(self.u), n))*dS(7+self.elec_list.index(electrode))
            # average is used since +/- sides of facet are arbitrary
            v = dolfin.assemble(m)
            cap_list[self.net_list.index(curr_net)] = cap_list[self.net_list.index(curr_net)] + v
        self.cap_matrix.append(cap_list)

    #get a delay estimate based on the simulation size and the capacitance
    def getDelays(self, bounds, temp):
        #have capacitance, resistivities, now need cross sectional area and length.
        #approximate length by characteristic length of simulation space (which was based on electrode placement)
        char_len = np.sqrt((bounds['xmax'] - bounds['xmin'])**2 + (bounds['ymax'] - bounds['ymin'])**2 +  (bounds['zmax'] - bounds['dielectric'])**2) # in angstrom
        print("char_len", char_len)
        #need the cross sectional area of the wire, maybe use minimum metal width ^2?
        mmw = 140 #angstroms
        print("mmw", mmw)
        rho = self.approxRes(temp) #in E-6 ohm cm
        rho *= 1E-6 #now in ohm cm
        rho *= 1E10/1E2 # now in ohm angstroms
        print("rho", rho)
        R = rho*char_len/mmw/mmw
        print("R", R)
        tau = R*max(max(self.cap_matrix))
        print("{:.2e}".format(1/tau/2/np.pi))

# print(res)
def test():
    rc = ResCap()
    print(rc.approxRes(77)) #Should be between 0.469 and 0.629

if __name__ == "__main__":
    test()
