 # @author:   Nathan
 # @created:  2019.01.28
 # @editted:  2019.01.29 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Plotting mechanisms for PoissonSolver

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from PIL import Image

class Plotter():
    #Constructor
    def __init__(self):
        pass

    #Produces a static potential figure with axes
    def saveAxesPotential(self,X,Y,Z,file_name):
        fig = plt.figure() #create a new figure
        plt.gca().invert_yaxis() # siqad has +x going right, +y going down
        max_val = np.max(np.abs(Z))
        # pin the colour scale extremes to +/- max value.
        norm = clrs.Normalize(vmin=-max_val, vmax=max_val)
        # set the colour map to use red/blue scheme
        plt.pcolormesh(X,Y,Z,norm=norm,cmap=plt.cm.get_cmap('RdBu_r'))
        cbar = plt.colorbar()
        cbar.set_label("Potential (V)")
        plt.xlabel("X (angstrom)") #set label text
        plt.ylabel("Y (angstrom)")
        #save the picture to file, and close
        plt.savefig(file_name, bbox_inches='tight', pad_inches=0)
        plt.close(fig)

    #Produces figures with the potential gradient, with axes
    def saveGrad(self,X,Y,Z,file_name):
        fig = plt.figure(frameon=False) #create a new figure
        plt.gca().invert_yaxis() # siqad has +x going right, +y going down
        Zgrad = np.gradient(Z) #get the gradient. Zgrad[0] is x direction, Zgrad[1] is y direction
        maxval = np.max(np.max(np.abs(Zgrad)))
        plt.title('Electric field from clocking electrodes')
        plt.quiver(X,Y,-Zgrad[1]/maxval,Zgrad[0]/maxval)
        plt.xlabel("X (angstrom)") #set label text
        plt.ylabel("Y (angstrom)")
        plt.savefig(file_name)
        plt.close(fig)

    #Produces a 2D slice of the potential at the given timestep without axes.
    def savePotential(self, X, Y, Z, step, file_name):
        fig = plt.figure(frameon=False)
        plt.gca().invert_yaxis()
        plt.axis('off')
        maxval = np.max(np.abs(Z))
        norm = clrs.Normalize(vmin=-maxval, vmax=maxval)
        plt.pcolormesh(X,Y,Z,norm=norm,cmap=plt.cm.get_cmap('RdBu_r'))
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.savefig(file_name)
        plt.close(fig)

    #Creates a gif from all "key{}" files in the given directory.
    def makeGif(self, mode, dir, key):
        dir_files = os.listdir(dir)
        if mode == "clock":
            images = []
            image_files = []
            for file in dir_files:
                if file.startswith(key):
                    image_files.append(os.path.join(dir, file))
            image_files.sort()
            for image_name in image_files:
                images.append(Image.open(image_name))
            images[0].save(os.path.join(dir, "SiAirBoundary.gif"),
                       save_all=True,
                       append_images=images[1:],
                       delay=0.5,
                       loop=0)

    def plotPotential(self, step, data, out_dir):
        X, Y, Z, nx, ny = data
        if step == 0:
            plot_file_name = os.path.join(out_dir,"SiAirPlot.png")
            self.saveAxesPotential(X, Y, Z, plot_file_name)
            grad_file_name = os.path.join(out_dir,'grad.pdf')
            self.saveGrad(X,Y,Z,grad_file_name)
        pot_file_name = os.path.join(out_dir,'SiAirBoundary{:03d}.png'.format(step))
        self.savePotential(X,Y,Z,step,pot_file_name)
