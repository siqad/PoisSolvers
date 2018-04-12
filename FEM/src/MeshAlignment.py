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
//PARAMETERS
mesh_resolution = 0.1;
x_max = 1.0;
x_min = 0.0;
y_max = 1.0;
y_min = 0.0;
z_max = 1.0;
z_min = 0.0;

// POINTS
Point(1) = {x_min, y_min, z_min, mesh_resolution}; //0, 0, 0
Point(2) = {x_min, y_min, z_max, mesh_resolution}; //0, 0, 1
Point(3) = {x_min, y_max, z_min, mesh_resolution}; //0, 1, 0
Point(4) = {x_min, y_max, z_max, mesh_resolution}; //0, 1, 1
Point(5) = {x_max, y_min, z_min, mesh_resolution}; //1, 0, 0
Point(6) = {x_max, y_min, z_max, mesh_resolution}; //1, 0, 1
Point(7) = {x_max, y_max, z_min, mesh_resolution}; //1, 1, 0
Point(8) = {x_max, y_max, z_max, mesh_resolution}; //1, 1, 1

//2D ITEMS
//LINES
//x = x_min surface
Line(1) = {1,2}; //0, 0, 0 -> 0, 0, 1
Line(2) = {2,4}; //0, 0, 1 -> 0, 1, 1
Line(3) = {4,3}; //0, 1, 1 -> 0, 1, 0
Line(4) = {3,1}; //0, 1, 0 -> 0, 0, 0

//x = x_max surface
Line(5) = {5,6}; //1, 0, 0 -> 1, 0, 1
Line(6) = {6,8}; //1, 0, 1 -> 1, 1, 1
Line(7) = {8,7}; //1, 1, 1 -> 1, 1, 0
Line(8) = {7,5}; //1, 1, 0 -> 1, 0, 0

//lines from x_max to x_min
Line(9) = {1,5}; //0, 0, 0 -> 1, 0, 0
Line(10) = {2,6}; //0, 0, 1 -> 1, 0, 1
Line(11) = {4,8}; //0, 1, 1 -> 1, 1, 1
Line(12) = {3,7}; //0, 1, 0 -> 1, 1, 0

//SURFACES
//x = x_min
Line Loop(13) = {1,2,3,4};
Plane Surface(14) = {13};
//x = x_max
Line Loop(15) = {5,6,7,8};
Plane Surface(16) = {15};
//y = y_min
Line Loop(17) = {9,5,-10,-1};
Plane Surface(18) = {17};
//y = y_max
Line Loop(19) = {12,-7,-11,3};
Plane Surface(20) = {19};
//z = z_min
Line Loop(21) = {9,-8,-12,4};
Plane Surface(22) = {21};
//z = z_max
Line Loop(23) = {10,6,-11,-2};
Plane Surface(24) = {23};

//CRACKS
Point(9) = {0.1, 0.1, 0.1, mesh_resolution};
Point(10) = {0.2, 0.2, 0.1, mesh_resolution}; 
Point(11) = {0.3, 0.1, 0.1, mesh_resolution};
Line(25) = {9,10};
Line(26) = {10,11};
Line(27) = {11,9};
Line Loop(28) = {25, 26, 27};
Plane Surface(29) = {28}
//Physical Surface(1) = {29};

//Volume
Surface Loop(30) = {14,16,18,20,22,24};
Volume(31) = {30};
Physical Volume(1) = {31};
Surface{29} In Volume{1}



"""

# Point(1) = {x_max,y_max,z_max,mesh_resolution};
# Point(2) = {x_max,y_max,z_min,mesh_resolution};
# Point(3) = {x_min,y_max,z_max,mesh_resolution};
# Point(4) = {x_min,y_min,z_max,mesh_resolution};
# Point(5) = {x_max,y_min,z_max,mesh_resolution};
# Point(6) = {x_max,y_min,z_min,mesh_resolution};
# Point(7) = {x_min,y_max,z_min,mesh_resolution};
# Point(8) = {x_min,y_min,z_min,mesh_resolution};
# Line(1) = {3,1};
# Line(2) = {3,7};
# Line(3) = {7,2};
# Line(4) = {2,1};
# Line(5) = {1,5};
# Line(6) = {5,4};
# Line(7) = {4,8};
# Line(8) = {8,6};
# Line(9) = {6,5};
# Line(10) = {6,2};
# Line(11) = {3,4};
# Line(12) = {8,7};
# Line Loop(13) = {-6,-5,-1,11};
# Plane Surface(14) = {13};
# Line Loop(15) = {4,5,-9,10};
# Plane Surface(16) = {15};
# Line Loop(17) = {-3,-12,8,10};
# Plane Surface(18) = {17};
# Line Loop(19) = {7,12,-2,11};
# Plane Surface(20) = {19};
# Line Loop(21) = {-4,-3,-2,1};
# Plane Surface(22) = {21};
# Line Loop(23) = {8,9,6,7};
# Plane Surface(24) = {23};
# 
# Surface Loop(25) = {14,24,-18,22,16,-20};
# Volume(26) = {25};
# Physical Volume(0) = {26};
# 
# """
with open('/home/nathan/git/PoisSolvers/FEM/data/domain.geo', 'w') as f: f.write(domain)

subprocess.call(['gmsh -3 /home/nathan/git/PoisSolvers/FEM/data/domain.geo'], shell=True)
subprocess.call(['dolfin-convert /home/nathan/git/PoisSolvers/FEM/data/domain.msh domain.xml'], shell=True)

mesh = Mesh('domain.xml')
# facet_f = MeshFunction('size_t', mesh, 'domain_facet_region.xml')

# for tag in (1, 2, 3):
    # print sum(1 for _ in SubsetIterator(facet_f, tag)), 'edges marked as', tag
plt.figure()
plot(mesh)
plt.show()
# interactive()
