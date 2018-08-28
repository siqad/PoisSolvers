 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Class definition for physics engine

import mesh_writer_3D as mw
import numpy as np
import os
import subdomains as sd
import dolfin

class PoissonSolver():
    def __init__(self, bounds):
        self.bounds = bounds
        self.mw = mw.MeshWriter()
        self.fields = []

    def setPaths(self, in_path="", out_path=""):
        if in_path != "":
            self.in_path = in_path
            self.abs_in_dir = os.path.abspath(os.path.dirname(self.in_path))
        if out_path != "":
            self.out_path = out_path
            self.abs_out_dir = os.path.abspath(os.path.dirname(self.out_path))

    def setResolution(self, res_scale):
        self.mw.resolution = min((self.bounds[1]-self.bounds[0]), \
                                 (self.bounds[3]-self.bounds[2]), \
                                 (self.bounds[5]-self.bounds[4])) / res_scale

    def createOuterBounds(self, resolution):
        self.mw.addBox([self.bounds[0],self.bounds[2],self.bounds[4]], \
                       [self.bounds[1],self.bounds[3],self.bounds[5]], resolution, "bound")

    def addDielectricSurface(self, resolution):
        self.mw.addSurface([self.bounds[0]+0.01*np.abs(self.bounds[0]),self.bounds[2]+0.01*np.abs(self.bounds[2]),self.bounds[6]],\
                      [self.bounds[1]-0.01*np.abs(self.bounds[1]),self.bounds[2]+0.01*np.abs(self.bounds[2]),self.bounds[6]],\
                      [self.bounds[1]-0.01*np.abs(self.bounds[1]),self.bounds[3]-0.01*np.abs(self.bounds[3]),self.bounds[6]],\
                      [self.bounds[0]+0.01*np.abs(self.bounds[0]),self.bounds[3]-0.01*np.abs(self.bounds[3]),self.bounds[6]],\
                      resolution, "seam")

    def addDielectricField(self, res_in, res_out):
        self.fields = []
        self.fields += [self.mw.addBoxField(res_in, res_out, \
                  [self.bounds[0]-0.001*np.abs(self.bounds[0]), \
                   self.bounds[1]+0.001*np.abs(self.bounds[1])], \
                  [self.bounds[2]-0.001*np.abs(self.bounds[2]), \
                   self.bounds[3]+0.001*np.abs(self.bounds[3])], \
                  [self.bounds[6]-0.05*np.abs(self.bounds[5]), \
                   self.bounds[6]+0.05*np.abs(self.bounds[5])])]
        self.fields = [self.mw.addMinField(self.fields)]

    def addElectrode(self, xs, ys, zs, resolution):
        self.mw.addBox([xs[0],ys[0],zs[0]], \
                  [xs[1],ys[1],zs[1]], resolution, "seam")

        #make resolution inside electrodes coarse
        self.fields += [self.mw.addBoxField(1.0, 0.0, \
                  [xs[0], xs[1]], [ys[0], ys[1]], [zs[0], zs[1]])]
        self.fields = [self.mw.addMaxField(self.fields)]
        self.fields += [self.mw.addBoxField(0.1, 1.0, \
                  [2.0*xs[0], 2.0*xs[1]], \
                  [2.0*ys[0], 2.0*ys[1]], \
                  [2.0*zs[0], 2.0*(zs[1])])]
        self.fields = [self.mw.addMinField(self.fields)]

    def addElectrodePoly(self, vertices, zs, resolution):
        self.mw.addPolygonVolume(vertices, zs, resolution)

    def setBGField(self, delta):
        bg_field_ind = self.mw.addMeanField(self.fields, delta)
        self.mw.setBGField(bg_field_ind)

    def writeGeoFile(self):
        # abs_in_dir = os.path.abspath(os.path.dirname(self.in_path))
        with open(os.path.join(self.abs_in_dir, 'domain.geo'), 'w') as f: f.write(self.mw.file_string)

    def setSubdomains(self):
        self.left = sd.Left(self.bounds[0]) #x
        self.top = sd.Top(self.bounds[3]) #y
        self.right = sd.Right(self.bounds[1]) #x
        self.bottom = sd.Bottom(self.bounds[2]) #y
        self.front = sd.Front(self.bounds[5]) #z
        self.back = sd.Back(self.bounds[4]) #z
        self.air = sd.Air((self.bounds[6], self.bounds[5]))

    def setElectrodeSubdomains(self, elec_list, zs):
        self.electrode = []
        for i in range(len(elec_list)):
            self.electrode.append(sd.Electrode([elec_list[i].x1, elec_list[i].x2], \
                                          [elec_list[i].y1, elec_list[i].y2], \
                                          zs ) )
            self.addElectrode([elec_list[i].x1, elec_list[i].x2], \
                            [elec_list[i].y1, elec_list[i].y2], \
                            zs, resolution=1.0)

    def setElectrodePolySubdomains(self, elec_poly_list, zs):
        self.electrode_poly = []
        for elec_poly in elec_poly_list:
            self.electrode_poly.append(sd.ElectrodePoly(elec_poly.vertex_list, \
                zs))
            ps.addElectrodePoly(elec_poly.vertex_list, zs, resolution=0.1)

    def markDomains(self, mesh):
        # Initialize mesh function for interior domains
        self.domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
        self.domains.set_all(0)
        self.air.mark(self.domains, 1)

    def markBoundaries(self, mesh, elec_list, elec_poly_list):
        self.boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
        self.boundaries.set_all(0)
        self.left.mark(self.boundaries, 1)
        self.top.mark(self.boundaries, 2)
        self.right.mark(self.boundaries, 3)
        self.bottom.mark(self.boundaries, 4)
        self.front.mark(self.boundaries, 5)
        self.back.mark(self.boundaries, 6)
        for i in range(len(elec_list)):
            self.electrode[i].mark(self.boundaries, 7+i)
        for i in range(len(elec_poly_list)):
            self.electrode_poly[i].mark(self.boundaries, 7+len(elec_list)+i)
