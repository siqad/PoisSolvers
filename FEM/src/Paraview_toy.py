import dolfin
import matplotlib.pyplot as plt
import mshr
from paraview import simple as pvs

box = mshr.Box(dolfin.Point(0,0,0), dolfin.Point(1,1,1))
sphere = mshr.Sphere(dolfin.Point(0.0,0.5,0.5), 0.25, 50)

combined = box + sphere

mesh = mshr.generate_mesh(combined, 50)

print mesh.topology().dim()

dolfin.File("/home/nathan/git/PoisSolvers/FEM/data/u.pvd") << mesh

wireframe = pvs.PVDReader(FileName="/home/nathan/git/PoisSolvers/FEM/data/u.pvd")
pvs.Show(wireframe)
prop = pvs.GetDisplayProperties(wireframe)
# print prop.GetProperty("Representation").Available
prop.Representation = "Wireframe"
# prop.Representation = "Surface with edges"

pvs.Interact(view=None)
pvs.Render()
# plt.figure()
# dolfin.plot(mesh, wireframe=True)
# plt.show()
