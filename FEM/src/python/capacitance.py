import numpy as np
import dolfin
import itertools

class CapacitanceEstimator():
    def __init__(self, t = None, r = None):
        self.cap_matrix = []
        self.dir = ""

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
        # print(cap_matrix)
        return cap_matrix

    def calcCaps(self):
        x0, x1, x2 = dolfin.MeshCoordinates(self.mesh)
        eps = dolfin.conditional(x2 <= 0.0, self.EPS_SI, self.EPS_DIELECTRIC)
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
