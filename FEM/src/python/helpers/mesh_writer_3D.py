 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Functions that create the GMSH .geo file

import matplotlib.pyplot as plt
import subprocess

class MeshWriter():
    def __init__(self, resolution=0.1):
        self.resolution = resolution
        self.file_string = ""
        self.ind_point = 0
        self.ind_2d = 0
        self.ind_bounding_vol = 0
        self.ind_vol = 0
        self.ind_phys_vol = 0
        self.ind_phys = 1
        self.ind_field = 1
        self.ind_poly = 0
        self.ind_boundaries = [[]]

    def addPolygonVolume(self, vertices, zs, scale):
        num_points = len(vertices)
        for vertex in vertices:
            self.addPoint([vertex[0],vertex[1],zs[0]], scale)
        line_inds = []
        for i in range(len(vertices)):
            if i == (len(vertices)-1) :
                line = self.addLineByIndex(self.ind_point-1-i, self.ind_point-1, scale)
            else:
                line = self.addLineByIndex(self.ind_point-1-i, self.ind_point-2-i, scale)
            line_inds.append(line)

        self.file_string += "Line Loop(%d) = {"%(self.ind_2d)
        for ind in line_inds:
            self.file_string += "%d"%(ind)
            if ind != line_inds[-1]:
                self.file_string += ","
        self.file_string += "};\n";
        self.ind_2d += 1

        self.file_string += "Plane Surface(%d) = {%d};\n"\
            %(self.ind_2d, self.ind_2d-1)
        self.ind_2d += 1
        self.file_string += "poly%d[] = Extrude {0,0,%.15f} { Surface{%d};};\n"\
            %(self.ind_poly, zs[1], self.ind_2d-1)
        self.file_string += "Surface{%d} In Volume{%d};\n"\
            %(self.ind_2d-1, self.ind_bounding_vol)
        walls = range(len(vertices)+2)
        #increase for extra surfaces from extrude
        self.ind_2d += len(vertices) + 1
        #increase for extra lines from extrude
        self.ind_2d += 5*len(vertices)
        self.ind_point += 4*len(vertices)
        for wall in walls:
            if wall != 1:
                self.file_string += "Surface{poly%d[%d]} In Volume{%d};\n"\
                    %(self.ind_poly, wall, self.ind_bounding_vol)
        self.file_string += "Delete { Volume{poly%d[1]}; }\n"%(self.ind_poly)
        self.ind_poly += 1

    #add a point at p, p in form [x, y]
    def addPoint(self, p, scale):
        self.file_string += "Point(%d) = {%.15f, %.15f, %.15f, %.15f};\n"\
            %(self.ind_point, p[0], p[1], p[2], self.resolution*scale)
        self.ind_point += 1
        return self.ind_point -1

    #line goes from p1->p2.
    def addLine(self, p1, p2, scale):
        self.addPoint(p1, scale)
        self.addPoint(p2, scale)
        self.file_string += "Line(%d) = {%d, %d};\n"\
            %(self.ind_2d, self.ind_point-2, self.ind_point-1)
        self.ind_2d += 1
        return self.ind_2d -1
    #add a line using existing points
    def addLineByIndex(self, p1_ind, p2_ind, scale):
        self.file_string += "Line(%d) = {%d, %d};\n"\
            %(self.ind_2d, p1_ind, p2_ind)
        self.ind_2d += 1
        return self.ind_2d -1


    def addLineLoop(self, a, b, c, d, option="bound", box=False):
        self.file_string += "Line Loop(%d) = {%d,%d,%d,%d};\n"\
            %(self.ind_2d,a,b,c,d)
        self.ind_2d += 1
        self.file_string += "Plane Surface(%d) = {%d};\n"\
            %(self.ind_2d,self.ind_2d-1)
        if option == "seam":
            self.file_string += "Surface{%d} In Volume{%d};\n"\
                %(self.ind_2d,self.ind_bounding_vol)
            if box == True:
                self.ind_boundaries[-1].append(self.ind_2d)
        if option == "bound":
            self.file_string += "Physical Surface(%d) = {%d};\n"\
                %(self.ind_phys, self.ind_2d)
            self.ind_phys += 1
        self.ind_2d += 1
        return self.ind_2d - 1
    #creates the outer boundary of the mesh
    #p in form [x,y]
    #start from bottom left, go counterclockwise
    #not enough that lines defined by points at same location, need to have SAME shared point index.
    def addBox(self, p1, p2, scale, option="bound"):
        x_min = min(p1[0],p2[0])
        y_min = min(p1[1],p2[1])
        z_min = min(p1[2],p2[2])
        x_max = max(p1[0],p2[0])
        y_max = max(p1[1],p2[1])
        z_max = max(p1[2],p2[2])

        self.addPoint([x_min,y_min,z_min], scale)
        self.addPoint([x_max,y_min,z_min], scale)
        self.addPoint([x_max,y_max,z_min], scale)
        self.addPoint([x_min,y_max,z_min], scale)

        self.addPoint([x_min,y_min,z_max], scale)
        self.addPoint([x_max,y_min,z_max], scale)
        self.addPoint([x_max,y_max,z_max], scale)
        self.addPoint([x_min,y_max,z_max], scale)

        l12 = self.addLineByIndex(self.ind_point-8, self.ind_point-7, scale)
        l11 = self.addLineByIndex(self.ind_point-7, self.ind_point-6, scale)
        l10 = self.addLineByIndex(self.ind_point-6, self.ind_point-5, scale)
        l9 = self.addLineByIndex(self.ind_point-5, self.ind_point-8, scale)

        l8 = self.addLineByIndex(self.ind_point-4, self.ind_point-3, scale)
        l7 = self.addLineByIndex(self.ind_point-3, self.ind_point-2, scale)
        l6 = self.addLineByIndex(self.ind_point-2, self.ind_point-1, scale)
        l5 = self.addLineByIndex(self.ind_point-1, self.ind_point-4, scale)

        l4 = self.addLineByIndex(self.ind_point-4, self.ind_point-8, scale)
        l3 = self.addLineByIndex(self.ind_point-3, self.ind_point-7, scale)
        l2 = self.addLineByIndex(self.ind_point-2, self.ind_point-6, scale)
        l1 = self.addLineByIndex(self.ind_point-1, self.ind_point-5, scale)

        #z=min
        self.addLineLoop(l12,l11,l10,l9,option,True)
        #z=max
        self.addLineLoop(l8,l7,l6,l5,option,True)
        #x=min
        self.addLineLoop(l4,-l9,-l1,l5,option,True)
        #x=max
        self.addLineLoop(l3,l11,-l2,-l7,option,True)
        #y=min
        self.addLineLoop(l4,l12,-l3,-l8,option,True)
        #y=max
        self.addLineLoop(l1,-l10,-l2,l6,option,True)
        self.ind_boundaries.append([])

        if option == "bound":
            #Create the surface loop from the surfaces made, and define the volume
            self.file_string += "Surface Loop(%d) = {%d,%d,%d,%d,%d,%d};\n"\
                %(self.ind_2d,self.ind_2d-1,self.ind_2d-3,self.ind_2d-5,self.ind_2d-7,self.ind_2d-9,self.ind_2d-11)
            self.ind_2d += 1
            self.file_string += "Volume(%d) = {%d};\n"\
                %(self.ind_2d,self.ind_2d-1)
            self.ind_bounding_vol = self.ind_2d
            self.ind_2d += 1
            self.file_string += "Physical Volume(%d) = {%d};\n"\
                %(self.ind_phys_vol, self.ind_2d-1)
        return

    def addPhysicalVolume(self):
        self.file_string += "Physical Volume(%d) = {%d};\n"\
            %(self.ind_phys_vol, self.ind_bounding_vol)

    def addSurface(self, p1, p2, p3, p4, scale, option):
        self.addPoint(p1, scale)
        self.addPoint(p2, scale)
        self.addPoint(p3, scale)
        self.addPoint(p4, scale)
        l4 = self.addLineByIndex(self.ind_point-4, self.ind_point-3, scale)
        l3 = self.addLineByIndex(self.ind_point-3, self.ind_point-2, scale)
        l2 = self.addLineByIndex(self.ind_point-2, self.ind_point-1, scale)
        l1 = self.addLineByIndex(self.ind_point-1, self.ind_point-4, scale)

        self.addLineLoop(l4,l3,l2,l1,option)

    #creates the outer boundary, and sets it as a volume.
    def addOuterBound(self, p1, p2, scale, option):
        self.addBox(p1, p2, scale, option)
        return

    #returns index of threshold field
    #The threshold field sets a resolution based on the distance away from a line
    #From 0->dist_min, the resolution is resolution*res_min_scale
    #From dist_max->inf, the resolution is resolution*res_max_scale
    #In between, the resolution is linearly interpolated.
    def addTHField(self, res_min_scale, res_max_scale, dist_min, dist_max):
        self.file_string += "Field[%d] = Attractor;\n"%(self.ind_field)
        self.file_string += "Field[%d].NNodesByEdge = 1000;\n"%(self.ind_field)
        self.file_string += "Field[%d].FacesList = {%d};\n"%(self.ind_field, self.ind_2d-1)
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
        self.ind_field += 1
        return self.ind_field-1

    #Add a minimum field, based on the given list of fields. The minimum field applies the minimum resolution of all fields.
    def addMaxField(self, field_list):
        #create comma separated string of field indexes
        fields_string = ",".join(map(str, field_list))
        self.file_string += 'Field[%d] = Max;\n'%(self.ind_field)
        self.file_string += 'Field[%d].FieldsList = {%s};\n'%(self.ind_field, fields_string)
        self.ind_field += 1
        return self.ind_field-1

    #Add a minimum field, based on the given list of fields. The minimum field applies the minimum resolution of all fields.
    def addMeanField(self, field_list, delta):
        #create comma separated string of field indexes
        fields_string = ",".join(map(str, field_list))
        self.file_string += 'Field[%d] = Mean;\n'%(self.ind_field)
        self.file_string += 'Field[%d].IField = %s;\n'%(self.ind_field, fields_string)
        self.file_string += 'Field[%d].Delta = %.15f;\n'%(self.ind_field, delta)
        self.ind_field += 1
        return self.ind_field-1

    def setBGField(self, ind):
        self.file_string += 'Background Field = %d;\n'%(ind)

    def addBoxField(self, res_in_scale, res_out_scale, xs, ys, zs):
        self.file_string += "Field[%d] = Box;\n"%(self.ind_field)
        self.file_string += "Field[%d].VIn = %.15f;\n"%(self.ind_field, self.resolution*res_in_scale)
        self.file_string += "Field[%d].VOut = %.15f;\n"%(self.ind_field, self.resolution*res_out_scale)
        self.file_string += "Field[%d].XMin = %.15f;\n"%(self.ind_field, xs[0])
        self.file_string += "Field[%d].XMax = %.15f;\n"%(self.ind_field, xs[1])
        self.file_string += "Field[%d].YMin = %.15f;\n"%(self.ind_field, ys[0])
        self.file_string += "Field[%d].YMax = %.15f;\n"%(self.ind_field, ys[1])
        self.file_string += "Field[%d].ZMin = %.15f;\n"%(self.ind_field, zs[0])
        self.file_string += "Field[%d].ZMax = %.15f;\n"%(self.ind_field, zs[1])
        self.ind_field += 1
        return self.ind_field-1

    def finalize(self):
        self.file_string += 'Coherence;'
