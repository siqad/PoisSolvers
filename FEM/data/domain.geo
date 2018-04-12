
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
Point(9) = {0.1, 0.1, z_min, mesh_resolution};
Point(10) = {x_max, y_max, z_min, mesh_resolution}; 
Line(25) = {9, 10};
Physical Line(1) = {25};
Line{25} In Surface{22};

//Volume
Surface Loop(26) = {14,16,18,20,22,24};
Volume(27) = {26};
Physical Volume(2) = {27};
