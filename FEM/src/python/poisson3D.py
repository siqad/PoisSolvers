 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2017.08.23 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Main function for physics engine

import sys
import dolfin
import subprocess
import numpy as np
import mesh_writer_3D as mw
import poisson_class
import helpers
import time
import os
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import siqadconn
import subdomains as sd
from dolfin_utils.meshconvert import meshconvert

######################PREAMBLE
#Get the I/O locations
in_path = sys.argv[1]
out_path = sys.argv[2]

#Instatiate an empty PoissonSolver object
ps = poisson_class.PoissonSolver()

#Create a SiQADConnector object for I/O between the tool and the simulation engine
sqconn = siqadconn.SiQADConnector("PoisSolver", in_path, out_path)

#...and give it to PoissonSolver
ps.setConnector(sqconn)

#Initialize a bunch of stuff in PS now that we have the connector,
#including simulation parameters, defining mesh geometry, and I/O paths.
#Users should modify existing sim_params after initialize(), but before createMesh()
#to ensure their modications are reflected.
ps.initialize()

#Kick-start the mesh definition process. Defines the boundaries, and inserts
#electrodes surfaces into the mesh, adding coarseness guides as well.
#Finishes by writing the mesh definition to a .geo file, and converting it to
#a dolfin compatible .xml format.
ps.createMesh()

#Marks the subdomains and boundaries as top, left, bottom, right, dielectric, silicon, metal, etc.
#Also sets constants (elementary charge, permittivity, etc.)
ps.setupSim()

#At this point, users can override the spatial charge density (defaults to 0 everywhere).
#Using a dolfin expression could work, for symbolic definition of the charge density.
# ps.f = dolfin.Constant("0.0")

for step in range(ps.steps):
    #Setup the FEM solver by declaring the function spaces, asserting the potential at electrode
    #surfaces, defining the variational form, and setting the initial guess and other parameters.
    #Also sets up solver method, preconditioner, and other hyper parameters.
    ps.setupDolfinSolver(step)

    #Solve the linear variational problem.
    ps.solve()

    #Export the data (potentials at DB locations, capacitances, 2D slices, etc.).
    ps.export(step)

ps.makeGif()
ps.getCaps()
