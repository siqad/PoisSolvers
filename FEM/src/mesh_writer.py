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
        self.ind_field = 1

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

    #returns index of threshold field
    #The threshold field sets a resolution based on the distance away from a line
    #From 0->dist_min, the resolution is resolution*res_min_scale
    #From dist_max->inf, the resolution is resolution*res_max_scale
    #In between, the resolution is linearly interpolated.
    def addTHField(self, res_min_scale, res_max_scale, dist_min, dist_max):
        self.file_string += "Field[%d] = Attractor;\n"%(self.ind_field)
        self.file_string += "Field[%d].NNodesByEdge = 100;\n"%(self.ind_field)
        self.file_string += "Field[%d].EdgesList = {%d};\n"%(self.ind_field, self.ind_2d-1)
        self.ind_field += 1
        self.file_string += 'Field[%d] = Threshold;\n'%(self.ind_field)
        self.file_string += "Field[%d].IField = %d;\n"%(self.ind_field, self.ind_field-1)
        self.file_string += 'Field[%d].LcMin = '%(self.ind_field)+str(self.resolution*res_min_scale)+';\n'
        self.file_string += 'Field[%d].LcMax = '%(self.ind_field)+str(self.resolution*res_max_scale)+';\n'
        self.file_string += 'Field[%d].DistMin = %.15f;\n'%(self.ind_field, dist_min)
        self.file_string += 'Field[%d].DistMax = %.15f;\n'%(self.ind_field, dist_max)
        self.ind_field +=1
        return self.ind_field-1
    
    #Add a minimum field, based on the given list of fields. The minimum field applies the minimum resolution of all fields.
    def addMinField(self, field_list):
        #create comma separated string of field indexes
        fields_string = ",".join(map(str, field_list))
        self.file_string += 'Field[%d] = Min;\n'%(self.ind_field)
        self.file_string += 'Field[%d].FieldsList = {%s};\n'%(self.ind_field, fields_string)
        self.file_string += 'Background Field = %d;\n'%(self.ind_field)
        self.ind_field += 1

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