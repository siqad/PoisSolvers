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
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from PIL import Image

class PoissonSolver():
    def __init__(self, bounds):
        keys = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'dielectric']
        self.bounds = dict(zip(keys, bounds))
        # self.bounds = bounds
        self.mw = mw.MeshWriter()
        self.fields = []
        self.sim_params = None
        self.in_path = None
        self.out_path_path = None
        self.elec_list = None
        self.elec_poly_list = None
        self.metal_offset = None
        self.metal_thickness = None
        self.cap_matrix = []
        self.net_list = []

    def setMetals(self, m_off, m_thick):
        self.metal_offset = m_off
        self.metal_thickness = m_thick

    def setPaths(self, in_path="", out_path=""):
        if in_path != "":
            self.in_path = in_path
            self.abs_in_dir = os.path.abspath(os.path.dirname(self.in_path))
        if out_path != "":
            self.out_path = out_path
            self.abs_out_dir = os.path.abspath(os.path.dirname(self.out_path))

    def setSimParams(self, sim_params=None):
        if sim_params:
            self.sim_params = sim_params

    def setResolution(self):
        res_scale = float(self.sim_params["sim_resolution"])
        self.mw.resolution = min((self.bounds['xmax']-self.bounds['xmin']), \
                                 (self.bounds['ymax']-self.bounds['ymin']), \
                                 (self.bounds['zmax']-self.bounds['zmin'])) / res_scale

    def createOuterBounds(self, resolution):
        self.mw.addBox([self.bounds['xmin'],self.bounds['ymin'],self.bounds['zmin']], \
                       [self.bounds['xmax'],self.bounds['ymax'],self.bounds['zmax']], resolution, "bound")

    def addDielectricSurface(self, resolution):
        self.mw.addSurface([self.bounds['xmin']+0.01*np.abs(self.bounds['xmin']),self.bounds['ymin']+0.01*np.abs(self.bounds['ymin']),self.bounds['dielectric']],\
                      [self.bounds['xmax']-0.01*np.abs(self.bounds['xmax']),self.bounds['ymin']+0.01*np.abs(self.bounds['ymin']),self.bounds['dielectric']],\
                      [self.bounds['xmax']-0.01*np.abs(self.bounds['xmax']),self.bounds['ymax']-0.01*np.abs(self.bounds['ymax']),self.bounds['dielectric']],\
                      [self.bounds['xmin']+0.01*np.abs(self.bounds['xmin']),self.bounds['ymax']-0.01*np.abs(self.bounds['ymax']),self.bounds['dielectric']],\
                      resolution, "seam")

    def addDielectricField(self, res_in, res_out):
        self.fields = []
        self.fields += [self.mw.addBoxField(res_in, res_out, \
                  [self.bounds['xmin']-0.001*np.abs(self.bounds['xmin']), \
                   self.bounds['xmax']+0.001*np.abs(self.bounds['xmax'])], \
                  [self.bounds['ymin']-0.001*np.abs(self.bounds['ymin']), \
                   self.bounds['ymax']+0.001*np.abs(self.bounds['ymax'])], \
                  [self.bounds['dielectric']-0.05*np.abs(self.bounds['zmax']), \
                   self.bounds['dielectric']+0.05*np.abs(self.bounds['zmax'])])]
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
        self.left = sd.Left(self.bounds['xmin']) #x
        self.top = sd.Top(self.bounds['ymax']) #y
        self.right = sd.Right(self.bounds['xmax']) #x
        self.bottom = sd.Bottom(self.bounds['ymin']) #y
        self.front = sd.Front(self.bounds['zmax']) #z
        self.back = sd.Back(self.bounds['zmin']) #z
        self.air = sd.Air((self.bounds['dielectric'], self.bounds['zmax']))

    def setElectrodeSubdomains(self):
        self.electrode = []
        zs = [self.metal_offset, self.metal_offset+self.metal_thickness]
        for i in range(len(self.elec_list)):
            self.electrode.append(sd.Electrode([self.elec_list[i].x1, self.elec_list[i].x2], \
                                          [self.elec_list[i].y1, self.elec_list[i].y2], \
                                          zs ) )
            self.addElectrode([self.elec_list[i].x1, self.elec_list[i].x2], \
                            [self.elec_list[i].y1, self.elec_list[i].y2], \
                            zs, resolution=1.0)

    def setElectrodePolySubdomains(self):
        self.electrode_poly = []
        zs = [self.metal_offset, self.metal_offset+self.metal_thickness]
        for elec_poly in self.elec_poly_list:
            self.electrode_poly.append(sd.ElectrodePoly(elec_poly.vertex_list, \
                zs))
            self.addElectrodePoly(elec_poly.vertex_list, zs, resolution=0.1)

    def markDomains(self, mesh):
        # Initialize mesh function for interior domains
        self.domains = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
        self.domains.set_all(0)
        self.air.mark(self.domains, 1)

    def markBoundaries(self, mesh):
        self.boundaries = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim()-1)
        self.boundaries.set_all(0)
        self.left.mark(self.boundaries, 1)
        self.top.mark(self.boundaries, 2)
        self.right.mark(self.boundaries, 3)
        self.bottom.mark(self.boundaries, 4)
        self.front.mark(self.boundaries, 5)
        self.back.mark(self.boundaries, 6)
        for i in range(len(self.elec_list)):
            self.electrode[i].mark(self.boundaries, 7+i)
        for i in range(len(self.elec_poly_list)):
            self.electrode_poly[i].mark(self.boundaries, 7+len(self.elec_list)+i)

    def getElecPotential(self, elec_list, step, steps, i):
        chi_si = 4.05 #eV
        phi_gold = 5.1 #eV
        phi_bi = phi_gold - chi_si
        elec_str = "Electrode "+str(i)+" is "
        if elec_list[i].electrode_type == 1:
            elec_str += "clocked, "
            tot_phase = elec_list[i].phase + step*360/steps
            potential_to_set = elec_list[i].potential*np.sin( np.deg2rad(tot_phase) )
        else:
            elec_str += "fixed, "
            potential_to_set = elec_list[i].potential
        if self.metal_offset > self.bounds['dielectric']:
            elec_str += "and above the dielectric interface."
        else:
            potential_to_set += phi_bi
            elec_str += "and below the dielectric interface."
        print(elec_str)
        return potential_to_set

    def getBoundaryComponent(self, u, v, ds):
        if self.sim_params["bcs"] == "robin":
            h_L = dolfin.Constant("0.0")
            h_R = dolfin.Constant("0.0")
            h_T = dolfin.Constant("0.0")
            h_Bo = dolfin.Constant("0.0")
            h_F = dolfin.Constant("0.0")
            h_Ba = dolfin.Constant("0.0")
            component =  h_L*u*v*ds(1) + h_R*u*v*ds(3) \
                + h_T*u*v*ds(2) + h_Bo*u*v*ds(4) \
                + h_F*u*v*ds(5) + h_Ba*u*v*ds(6)
        elif self.sim_params["bcs"] == "neumann":
            g_L = dolfin.Constant("0.0")
            g_R = dolfin.Constant("0.0")
            g_T = dolfin.Constant("0.0")
            g_Bo = dolfin.Constant("0.0")
            g_F = dolfin.Constant("0.0")
            g_Ba = dolfin.Constant("0.0")
            component = - g_L*v*ds(1) - g_R*v*ds(3) \
                 - g_T*v*ds(2) - g_Bo*v*ds(4) \
                 - g_F*v*ds(5) - g_Ba*v*ds(6)
        return component

    def getSteps(self):
        mode = str(self.sim_params["mode"])
        if mode == "standard":
            steps = 1
        elif mode == "clock":
            steps = int(self.sim_params["steps"])
        elif mode == "cap":
            steps = len(self.net_list)
        return steps

    def createNetlist(self):
        for i in range(len(self.elec_list)):
            if self.elec_list[i].net not in self.net_list:
                self.net_list.append(self.elec_list[i].net)
        for i in range(len(self.elec_poly_list)):
            if self.elec_poly_list[i].net not in self.net_list:
                self.net_list.append(self.elec_poly_list[i].net)

    def setElectrodePotentials(self, step, steps, V):
        mode = str(self.sim_params["mode"])
        if mode == "standard" or mode == "clock":
            for i in range(len(self.elec_list)):
                potential_to_set = self.getElecPotential(self.elec_list, step, steps, i)
                self.bcs.append(dolfin.DirichletBC(V, float(potential_to_set), self.boundaries, 7+i))
            for i in range(len(self.elec_poly_list)):
                potential_to_set = self.getElecPotential(self.elec_poly_list, step, steps, i)
                self.bcs.append(dolfin.DirichletBC(V, float(potential_to_set), self.boundaries, 7+len(self.elec_list)+i))
        elif mode == "cap":
            # for i in range(len(self.elec_list)+len(self.elec_poly_list)):
            #     if i == step:
            #         self.bcs.append(dolfin.DirichletBC(V, float(1.0), self.boundaries, 7+i))
            #     else:
            #         self.bcs.append(dolfin.DirichletBC(V, float(0.0), self.boundaries, 7+i))
            for i in range(len(self.elec_list)):
                if self.net_list[step] == self.elec_list[i].net:
                    self.bcs.append(dolfin.DirichletBC(V, float(1.0), self.boundaries, 7+i))
                else:
                    self.bcs.append(dolfin.DirichletBC(V, float(0.0), self.boundaries, 7+i))
            for i in range(len(self.elec_poly_list)):
                if self.net_list[step] == self.elec_poly_list[i].net:
                    self.bcs.append(dolfin.DirichletBC(V, float(1.0), self.boundaries, 7+len(self.elec_list)+i))
                else:
                    self.bcs.append(dolfin.DirichletBC(V, float(0.0), self.boundaries, 7+len(self.elec_list)+i))

    def saveAxesPotential(self,X,Y,Z,filename):
        fig = plt.figure()
        plt.gca().invert_yaxis()
        maxval = np.max(np.abs(Z))
        norm = clrs.Normalize(vmin=-maxval, vmax=maxval)
        plt.pcolormesh(X,Y,Z,norm=norm,cmap=plt.cm.get_cmap('RdBu_r'))
        cbar = plt.colorbar()
        cbar.set_label("Potential (V)")
        locs, labels = plt.yticks()
        labels = []
        for loc in locs:
            labels += [str(round(loc*1e9, 2))]
        plt.yticks(locs, labels)
        locs, labels = plt.xticks()
        labels = []
        for loc in locs:
            labels += [str(round(loc*1e9, 2))]
        plt.xticks(locs, labels)
        plt.xlabel("X (nm)")
        plt.ylabel("Y (nm)")
        savestring = os.path.join(self.abs_out_dir,filename)
        plt.savefig(savestring, bbox_inces="tight", pad_inches=0)
        plt.close(fig)

    def saveGrad(self, X, Y, Z, index):
        fig = plt.figure(frameon=False)
        plt.gca().invert_yaxis()
        Zgrad = np.gradient(Z)
        maxval = np.max(np.abs(Zgrad[index]))
        norm = clrs.Normalize(vmin=-maxval, vmax=maxval)
        plt.pcolormesh(X,Y,Zgrad[index],norm=norm,cmap=plt.cm.get_cmap('RdBu_r'))
        cbar = plt.colorbar()
        cbar.set_label("E field (V/m)")
        savestring = os.path.join(self.abs_out_dir,'grad{}.png'.format(index))
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.savefig(savestring)
        plt.close(fig)

    def savePotential(self, X, Y, Z, step):
        fig = plt.figure(frameon=False)
        plt.gca().invert_yaxis()
        plt.axis('off')
        maxval = np.max(np.abs(Z))
        norm = clrs.Normalize(vmin=-maxval, vmax=maxval)
        plt.pcolormesh(X,Y,Z,norm=norm,cmap=plt.cm.get_cmap('RdBu_r'))
        savestring = os.path.join(self.abs_out_dir,'SiAirBoundary{:03d}.png'.format(step))
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.savefig(savestring)
        plt.close(fig)

    def create2DSlice(self, u):
        depth = float(self.sim_params['slice_depth'])*1e-9
        nx = int(self.sim_params['image_resolution'])
        ny = nx
        x = np.linspace(self.bounds['xmin'], self.bounds['xmax'], nx)
        y = np.linspace(self.bounds['ymin'], self.bounds['ymax'], ny)
        X, Y = np.meshgrid(x, y)
        z = np.array([u(i, j, self.bounds['dielectric']-depth) for j in y for i in x])
        Z = z.reshape(nx, ny)
        return X, Y, Z, nx, ny

    def makeGif(self):
        mode = str(self.sim_params["mode"])
        if mode == "clock":
            images = []
            image_files = []
            for file in os.listdir(os.path.dirname(self.in_path)):
                if file.startswith("SiAirBoundary"):
                    image_files.append(os.path.join(self.abs_in_dir, file))
            image_files.sort()
            for image_name in image_files:
                images.append(Image.open(image_name))
            images[0].save(os.path.join(self.abs_in_dir, "SiAirBoundary.gif"),
                       save_all=True,
                       append_images=images[1:],
                       delay=0.5,
                       loop=0)

    def getCaps(self):
        mode = str(self.sim_params["mode"])
        if mode == "cap":
            cap_matrix = np.array(self.cap_matrix)
            for i in range(len(cap_matrix)):
                tot_cap = 0
                for cap in cap_matrix[i]:
                    tot_cap = tot_cap+cap
                print("C_net{} = {}F".format(self.net_list[i],tot_cap))
            print(cap_matrix)

    def calcCaps(self, u, mesh, EPS_SI, EPS_AIR):
        mode = str(self.sim_params["mode"])
        if mode == "cap":
            x0, x1, x2 = dolfin.MeshCoordinates(mesh)
            eps = dolfin.conditional(x2 <= 0.0, EPS_SI, EPS_AIR)
            cap_list = [0.0]*len(self.net_list)
            for i in range(len(self.elec_list)):
                curr_net = self.elec_list[i].net
                dS = dolfin.Measure("dS")[self.boundaries]
                n = dolfin.FacetNormal(mesh)
                m = dolfin.avg(dolfin.dot(eps*dolfin.grad(u), n))*dS(7+i)
                # average is used since +/- sides of facet are arbitrary
                v = dolfin.assemble(m)
                print("\int grad(u) * n ds({}) = ".format(7+i), v)
                print(self.net_list.index(curr_net))
                cap_list[self.net_list.index(curr_net)] = cap_list[self.net_list.index(curr_net)] + v
                print(cap_list)
            for i in range(len(self.elec_poly_list)):
                curr_net = self.elec_poly_list[i].net
                dS = dolfin.Measure("dS")[self.boundaries]
                n = dolfin.FacetNormal(mesh)
                m = dolfin.avg(dolfin.dot(eps*dolfin.grad(u), n))*dS(7+len(self.elec_list)+i)
                # average is used since +/- sides of facet are arbitrary
                v = dolfin.assemble(m)
                print("\int grad(u) * n ds({}) = ".format(7+len(self.elec_list)+i), v)
                print(self.net_list.index(curr_net))
                cap_list[self.net_list.index(curr_net)] = cap_list[self.net_list.index(curr_net)] + v
                print(cap_list)
            self.cap_matrix.append(cap_list)
