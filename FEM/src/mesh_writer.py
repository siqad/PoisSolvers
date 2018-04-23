import matplotlib.pyplot as plt
import dolfin
import subprocess

class MeshWriter():
    def __init__(self, resolution=0.1):
        self.resolution = resolution
        self.file_string = ""
        self.ind_point = 0
        self.ind_2d = 0
        self.ind_surface = 0
        self.ind_phys = 0
        self.ind_phys_surface = 0

    #add a point at p, p in form [x, y]
    def addPoint(self, p, scale):
        self.file_string += "Point(%d) = {%.15f, %.15f, 0.000, %.15f};\n"\
            %(self.ind_point, p[0], p[1], self.resolution*scale)
        self.ind_point += 1
        return
    
    #add a point at p, p in form [x, y]
    def addPointToSurface(self, p, scale):
        self.file_string += "Point(%d) = {%.15f, %.15f, 0.000, %.15f};\n"\
            %(self.ind_point, p[0], p[1], self.resolution*scale)
        self.ind_point += 1
        self.file_string += "Point{%d} In Surface{%d};\n"\
            %(self.ind_point-1, self.ind_surface)
        return
    
    #line goes from p1->p2.
    def addLine(self, p1, p2, scale):
        self.addPoint(p1, scale)
        self.addPoint(p2, scale)
        self.file_string += "Line(%d) = {%d, %d};\n"\
            %(self.ind_2d, self.ind_point-2, self.ind_point-1)
        self.ind_2d += 1
        return
    #add a line using existing points
    def addLineByIndex(self, p1_ind, p2_ind, scale):
        self.file_string += "Line(%d) = {%d, %d};\n"\
            %(self.ind_2d, p1_ind, p2_ind)
        self.ind_2d += 1
        return
    
    #creates the outer boundary of the mesh
    #p in form [x,y]
    #start from bottom left, go counterclockwise
    #not enough that lines defined by points at same location, need to have SAME shared point index.
    def addRect(self, p1, p2, scale):
        x_min = min(p1[0],p2[0])
        y_min = min(p1[1],p2[1])
        x_max = max(p1[0],p2[0])
        y_max = max(p1[1],p2[1])
        self.addPoint([x_min,y_min], scale)
        self.addPoint([x_max,y_min], scale)
        self.addPoint([x_max,y_max], scale)
        self.addPoint([x_min,y_max], scale)
        self.addLineByIndex(self.ind_point-4, self.ind_point-3, scale)
        self.addLineByIndex(self.ind_point-3, self.ind_point-2, scale)
        self.addLineByIndex(self.ind_point-2, self.ind_point-1, scale)
        self.addLineByIndex(self.ind_point-1, self.ind_point-4, scale)
        self.file_string += "Line Loop(%d) = {%d,%d,%d,%d};\n"\
            %(self.ind_2d,self.ind_point-4,self.ind_point-3,self.ind_point-2,self.ind_point-1)
        self.ind_2d += 1
        self.file_string += "Plane Surface(%d) = {%d};\n"\
            %(self.ind_2d,self.ind_2d-1)
        self.ind_surface = self.ind_2d
        self.ind_2d += 1
        return

    #creates the outer boundary, and sets it as a surface.
    def addOuterBound(self, p1, p2, scale):
        self.addRect(p1, p2, scale)
        self.file_string += "Physical Surface(%d) = {%d};\n"\
            %(self.ind_phys_surface, self.ind_2d-1)
        return

    #adds a crack that the mesh must align to.
    def addCrack(self, p1, p2, scale):
        self.addLine(p1, p2, scale)
        self.file_string += "Physical Line(%d) = {%d};\n"\
            %(self.ind_phys, self.ind_2d-1)
        self.ind_phys += 1
        self.file_string += "Line{%d} In Surface{%d};\n"\
            %(self.ind_2d-1, self.ind_surface)        
        return
        
    def addCrackByIndex(self, p1_ind, p2_ind, scale):
        self.file_string += "Line(%d) = {%d, %d};\n"\
            %(self.ind_2d, p1_ind, p2_ind)
        self.ind_2d += 1
        self.file_string += "Physical Line(%d) = {%d};\n"\
            %(self.ind_phys,self.ind_2d-1)
        self.ind_phys += 1
        self.file_string += "Line{%d} In Surface{%d};\n"\
            %(self.ind_2d-1, self.ind_surface)        
        return
    
    #defines a box made from cracks, which the mesh will align itself to.
    #will set the internal resolution to standard,
    #keeping scaled resolution at the crack boundary
    def addCrackBox(self, p1, p2, scale):
        x_min = min(p1[0],p2[0])
        y_min = min(p1[1],p2[1])
        x_max = max(p1[0],p2[0])
        y_max = max(p1[1],p2[1])
        self.addPoint([x_min,y_min], scale)
        self.addPoint([x_max,y_min], scale)
        self.addPoint([x_max,y_max], scale)
        self.addPoint([x_min,y_max], scale)
        self.addCrackByIndex(self.ind_point-4, self.ind_point-3, scale)
        self.addCrackByIndex(self.ind_point-3, self.ind_point-2, scale)
        self.addCrackByIndex(self.ind_point-2, self.ind_point-1, scale)
        self.addCrackByIndex(self.ind_point-1, self.ind_point-4, scale)
        return
    
# mw = MeshWriter()
# 
# mw.addOuterBound([0.0,0.0], [1.0,1.0], 1)
# mw.addCrack([0.01,0.8],[0.99,0.8],0.1)
# mw.addCrack([0.01,0.7],[0.99,0.7],0.5)
# mw.addCrackBox([0.4,0.4],[0.6,0.6],0.1)
# # print mw.file_string
# 
# with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
# subprocess.call(['gmsh -2 ../data/domain.geo'], shell=True)
# subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
# mesh = dolfin.Mesh('domain.xml')
# plt.figure()
# dolfin.plot(mesh)
# plt.show()