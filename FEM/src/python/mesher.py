 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Class for meshing tools. Will handle the meshing portion of the solver.

import mesh_writer_3D as mw
import numpy as np
import os
import subdomains as sd
import subprocess

class Mesher():
    #Constructor
    def __init__(self):
        #mesh generation
        self.mw = mw.MeshWriter() #class for creating the mesh geometry
        self.resolution = None
        self.fields = [] #resolution fields for the mesh.
        self.bounds = None
        self.dir = None
        self.elec_list = []
        self.metal_params = []

    def createGeometry(self):
        self.setResolution()
        self.createOuterBounds(1.0)
        self.createSubstrateBox(1.0)
        #if we don't need the doping, then no need for extra resolution
        if self.eqn != "laplace":
            self.createDopingBox(0.25)
        surfaces = ['left', 'top', 'right', 'bottom', 'front', 'back', 'air']
        sd_list = list(self.setSubdomains())
        subdomains = dict(zip(surfaces, sd_list))
        electrodes = self.setElectrodeSubdomains()
        self.setBGField()
        self.finalize()
        self.writeGeoFile()
        self.runGMSH()
        return subdomains, electrodes


    def runGMSH(self, file_name="domain.geo"):
        subprocess.call(["gmsh", "-format", "msh22", "-3", os.path.join(self.dir,file_name)])

    def setResolution(self):
        x_min = self.bounds['xmin']
        y_min = self.bounds['ymin']
        z_min = self.bounds['zmin']
        x_max = self.bounds['xmax']
        y_max = self.bounds['ymax']
        z_max = self.bounds['zmax']
        # base the resolution on the average of the three dimensions
        self.mw.resolution = ((x_max-x_min) + (y_max-y_min) + (z_max-z_min)) / 3.0 / self.resolution

    def createOuterBounds(self, resolution=1.0):
        print("Creating outer boundaries...")
        self.mw.addBox([self.bounds['xmin'],self.bounds['ymin'],self.bounds['zmin']], \
                       [self.bounds['xmax'],self.bounds['ymax'],self.bounds['zmax']], self.mw.resolution*resolution, option="bound")

    def createSubstrateBox(self, resolution=1.0):
        #Add in a box on the negative half, to get the plane surface at z = 0
        self.mw.addBox([self.bounds['xmin'],self.bounds['ymin'],self.bounds['zmin']], \
                       [self.bounds['xmax'],self.bounds['ymax'],self.bounds['dielectric']], self.mw.resolution*resolution, option="seam")
        # z_size = self.bounds['dielectric'] - self.bounds['zmin']
        # self.fields += [self.mw.addBoxField(resolution, 1.0, \
        #           [self.bounds['xmin'], self.bounds['xmax']], \
        #           [self.bounds['ymin'], self.bounds['ymax']], \
        #           [self.bounds['zmin'], self.bounds['dielectric']])]
        # self.fields = [self.mw.addMinField(self.fields)]


    def createDopingBox(self, resolution=1.0):
        #Add in a box on the negative half, to get the plane surface at z = 0
        z_size = self.bounds['dielectric'] - self.bounds['zmin']
        self.mw.addBox([self.bounds['xmin'],self.bounds['ymin'],self.bounds['zmin']+0.25*z_size], \
                       [self.bounds['xmax'],self.bounds['ymax'],self.bounds['dielectric']-0.25*z_size], self.mw.resolution*resolution, option="seam")
        # z_size = self.bounds['dielectric'] - self.bounds['zmin']
        self.fields += [self.mw.addBoxField(resolution, 1.0, \
                  [self.bounds['xmin'], self.bounds['xmax']], \
                  [self.bounds['ymin'], self.bounds['ymax']], \
                  [self.bounds['zmin']+0.25*z_size, self.bounds['dielectric']-0.25*z_size])]
        self.fields = [self.mw.addMinField(self.fields)]

    def addElectrode(self, electrode, resolution):
        x_min = self.bounds['xmin']
        y_min = self.bounds['ymin']
        x_max = self.bounds['xmax']
        y_max = self.bounds['ymax']
        z_min = self.bounds['zmin']
        z_max = self.bounds['zmax']
        xs = [electrode.x1,electrode.x2]
        ys = [electrode.y1,electrode.y2]
        zs = [electrode.z1,electrode.z2]
        self.mw.addBox([xs[0],ys[0],zs[0]], \
                  [xs[1],ys[1],zs[1]], resolution,angle=electrode.angle,option="seam")
                  # [xs[1],ys[1],zs[1]], resolution,angle=electrode.angle,option="bound")
        #The physical extent of the field
        dist_x = 0.01*(x_max - x_min)
        dist_y = 0.01*(y_max - y_min)
        dist_z = 0.01*(z_max - z_min)

        self.fields += [self.mw.addBoxField(resolution, 1.0, \
                  [xs[0]-dist_x, xs[1]+dist_x], \
                  [ys[0]-dist_y, ys[1]+dist_y], \
                  [zs[0]-dist_z, zs[1]+dist_z])]
        self.fields = [self.mw.addMinField(self.fields)]

    def setBGField(self, delta=5):
        # bg_field_ind = self.mw.addMeanField(self.fields, delta)
        self.mw.setBGField(self.fields[0])

    def writeGeoFile(self):
        if self.dir == None:
            print("Directory not set for writing geometry file.")
        else:
            print("Writing mesh definition to file...")
            with open(os.path.join(self.dir, 'domain.geo'), 'w') as f: f.write(self.mw.file_string)

    def setSubdomains(self):
        print("Create subdomains and fields...")
        self.left = sd.Left(self.bounds['xmin'])
        self.top = sd.Top(self.bounds['ymax'])
        self.right = sd.Right(self.bounds['xmax'])
        self.bottom = sd.Bottom(self.bounds['ymin'])
        self.front = sd.Front(self.bounds['zmax'])
        self.back = sd.Back(self.bounds['zmin'])
        self.air = sd.Air((self.bounds['dielectric'], self.bounds['zmax']))
        return self.left, self.top, self.right, self.bottom, self.front, self.back, self.air

    def setElectrodeSubdomains(self):
        self.electrodes = []
        for elec in self.elec_list:
            layer_id = elec.layer_id # extract the layer id
            zs = [self.metal_params[layer_id][0], sum(self.metal_params[layer_id])]
            zs = [min(zs),max(zs)]
            elec.z1 = zs[0]
            elec.z2 = zs[1] #fill in the z dimension
            self.electrodes.append(sd.Electrode(elec)) #add to the list of electrodes
            self.addElectrode(elec, resolution=1.0) #add the electrode into the mesh
        return self.electrodes

    def finalize(self):
        self.mw.finalize();
