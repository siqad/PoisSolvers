import mesh_writer_3D as mw
import numpy as np

class PoissonSolver():
    def __init__(self, bounds):
        self.bounds = bounds
        self.mw = mw.MeshWriter()
        self.fields = []

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
