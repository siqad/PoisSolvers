import dolfin
import subprocess
import numpy as np
import time
import sys
sys.path.insert(0, "./python/")
import electrode_parser
import mesh_writer_3D as mw

print len(sys.argv)
print(sys.argv)
print "Hello World"

nx, ny = (10, 10)
x = np.linspace(2, 10, nx)
y = np.linspace(5, 20, ny)
yv, xv = np.meshgrid(y, x)
print "xv: ", xv
print "yv: ", yv