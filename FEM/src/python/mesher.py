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
        self.createOuterBounds()
        self.addDielectricSurface()
        self.addDielectricField()
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
        subprocess.call(["gmsh", "-3", os.path.join(self.dir,file_name)])

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
                       [self.bounds['xmax'],self.bounds['ymax'],self.bounds['zmax']], self.resolution, option="bound")

    def addDielectricSurface(self, resolution=1.0):
        print("Inserting dielectric surface...")
        x_min = self.bounds['xmin']
        y_min = self.bounds['ymin']
        x_max = self.bounds['xmax']
        y_max = self.bounds['ymax']
        z = self.bounds['dielectric']
        # add the surface using 4 points of a rectangle
        self.mw.addSurface([x_min,y_min,z], [x_max,y_min,z], [x_max,y_max,z], [x_min,y_max,z], 1, option="seam")

    def addDielectricField(self, res_in=0.25, res_out=1.0):
        #Shift the points in a bit, so the field doesn't affect the boundary vertices.
        x_length = abs(self.bounds['xmax']-self.bounds['xmin'])
        y_length = abs(self.bounds['ymax']-self.bounds['ymin'])
        # z_length = abs(self.bounds['zmax']-self.bounds['zmin'])
        x_min = self.bounds['xmin'] + 0.05*x_length
        y_min = self.bounds['ymin'] + 0.05*y_length
        x_max = self.bounds['xmax'] - 0.05*x_length
        y_max = self.bounds['ymax'] - 0.05*y_length
        dielec = self.bounds['dielectric']
        z_min = self.bounds['zmin']
        z_max = self.bounds['zmax']
        self.fields = []
        #construct a field using the given resolutions and the dimensions of the box.
        self.fields += [self.mw.addBoxField(res_in, res_out, [x_min, x_max], [y_min, y_max], \
                  [dielec-0.1*np.abs(z_min), dielec+0.1*np.abs(z_max)])]
                  # [dielec-0.1*z_length, dielec+0.1*z_length])]
        self.fields = [self.mw.addMinField(self.fields)]

    def addElectrode(self, electrode, resolution):
        xs = [electrode.x1,electrode.x2]
        ys = [electrode.y1,electrode.y2]
        zs = [electrode.z1,electrode.z2]
        self.mw.addBox([xs[0],ys[0],zs[0]], \
                  [xs[1],ys[1],zs[1]], resolution,angle=electrode.angle,option="seam")
        #The physical extent of the field
        dist_x = 0.5*(xs[1] - xs[0])
        dist_y = 0.5*(ys[1] - ys[0])
        dist_z = 0.5*(zs[1] - zs[0])

        self.fields += [self.mw.addBoxField(0.1, 1.0, \
                  [xs[0]-dist_x, xs[1]+dist_x], \
                  [ys[0]-dist_y, ys[1]+dist_y], \
                  [zs[0]-dist_z, zs[1]+dist_z])]
        self.fields = [self.mw.addMinField(self.fields)]

    def setBGField(self, delta=10):
        bg_field_ind = self.mw.addMeanField(self.fields, delta)
        self.mw.setBGField(bg_field_ind)

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
