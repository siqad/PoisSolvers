from dolfin import *
import subprocess
import matplotlib.pyplot as plt

# domain = '''
# //mesh_resolution: small->fine, large->coarse
# mesh_resolution = 0.05;
# Point(1) = {0, 0, 0, mesh_resolution};
# Point(2) = {1, 0, 0, mesh_resolution};
# Point(3) = {1, 1, 0, mesh_resolution};
# Point(4) = {0, 1, 0, mesh_resolution};
# // 2d domain
# Line(1) = {1, 2};
# Line(2) = {2, 3};
# Line(3) = {3, 4};
# Line(4) = {4, 1};
# Line Loop(5) = {1, 2, 3, 4};
# Plane Surface(6) = {5};
# 
# // Cracks are locations that the mesh is forced to align to
# 
# // 1
# Point(5) = {0.25, 0.25, 0, mesh_resolution};
# Point(6) = {0.5, 0.75, 0, mesh_resolution}; 
# Line(7) = {5, 6};
# Physical Line(1) = {7};
# Line{7} In Surface{6};
# 
# // 2
# Point(7) = {0.75, 0.25, 0, mesh_resolution};
# Line(8) = {6, 7};
# Physical Line(2) = {8};
# Line{8} In Surface{6};
# 
# Line(9) = {7, 5};
# Physical Line(3) = {9};
# Line{9} In Surface{6};
# 
# Physical Surface(1) = {6};
# '''

domain = """
lcar1 = .1;

length = 1.0;
height = 1.0;
depth = 1.0;

Point(newp) = {length/2,height/2,depth,lcar1}; /* Point      1 */
Point(newp) = {length/2,height/2,0,lcar1}; /* Point      2 */
Point(newp) = {-length/2,height/2,depth,lcar1}; /* Point      3 */
Point(newp) = {-length/2,-height/2,depth,lcar1}; /* Point      4 */
Point(newp) = {length/2,-height/2,depth,lcar1}; /* Point      5 */
Point(newp) = {length/2,-height/2,0,lcar1}; /* Point      6 */
Point(newp) = {-length/2,height/2,0,lcar1}; /* Point      7 */
Point(newp) = {-length/2,-height/2,0,lcar1}; /* Point      8 */
Line(1) = {3,1};
Line(2) = {3,7};
Line(3) = {7,2};
Line(4) = {2,1};
Line(5) = {1,5};
Line(6) = {5,4};
Line(7) = {4,8};
Line(8) = {8,6};
Line(9) = {6,5};
Line(10) = {6,2};
Line(11) = {3,4};
Line(12) = {8,7};
Line Loop(13) = {-6,-5,-1,11};
Plane Surface(14) = {13};
Line Loop(15) = {4,5,-9,10};
Plane Surface(16) = {15};
Line Loop(17) = {-3,-12,8,10};
Plane Surface(18) = {17};
Line Loop(19) = {7,12,-2,11};
Plane Surface(20) = {19};
Line Loop(21) = {-4,-3,-2,1};
Plane Surface(22) = {21};
Line Loop(23) = {8,9,6,7};
Plane Surface(24) = {23};

Surface Loop(25) = {14,24,-18,22,16,-20};
Volume(26) = {25};

"""
with open('/home/nathan/git/PoisSolvers/FEM/data/domain.geo', 'w') as f: f.write(domain)

subprocess.call(['gmsh -3 /home/nathan/git/PoisSolvers/FEM/data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert /home/nathan/git/PoisSolvers/FEM/data/domain.msh domain.xml'], shell=True)

mesh = Mesh('domain.xml')
facet_f = MeshFunction('size_t', mesh, 'domain_facet_region.xml')

for tag in (1, 2, 3):
    print sum(1 for _ in SubsetIterator(facet_f, tag)), 'edges marked as', tag
plt.figure()
plot(mesh)
plt.show()
# interactive()
