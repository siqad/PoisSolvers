import matplotlib.pyplot as plt
import dolfin
import subprocess

class MeshWriter():
    def __init__(self, resolution=0.1):
        self.resolution = resolution
        self.file_string = ""
        self.ind_point = 0
        self.ind_2d = 0
        self.ind_vol = 0
        self.ind_phys_vol = 0
        self.ind_phys = 1
        self.ind_field = 1
        self.ind_phys_surf = []

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


    def addLineLoop(self, a, b, c, d, option="bound"):
        self.file_string += "Line Loop(%d) = {%d,%d,%d,%d};\n"\
            %(self.ind_2d,a,b,c,d)
        self.ind_2d += 1
        self.file_string += "Plane Surface(%d) = {%d};\n"\
            %(self.ind_2d,self.ind_2d-1)
        if option == "seam":
            self.file_string += "Surface{%d} In Volume{%d};\n"\
                %(self.ind_2d,self.ind_vol)
        # self.ind_phys_surf.append(self.ind_2d)
        if option == "bound":
            self.file_string += "Physical Surface(%d) = {%d};\n"\
                %(self.ind_phys, self.ind_2d)
            self.ind_phys_vol += 1
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
        self.addLineLoop(l12,l11,l10,l9,option)
        #z=max
        self.addLineLoop(l8,l7,l6,l5,option)
        #x=min
        self.addLineLoop(l4,-l9,-l1,l5,option)
        #x=max
        self.addLineLoop(l3,l11,-l2,-l7,option)
        #y=min
        self.addLineLoop(l4,l12,-l3,-l8,option)
        #y=max
        self.addLineLoop(l1,-l10,-l2,l6,option)

        if option == "bound":
            #Create the surface loop from the surfaces made, and define the volume
            self.file_string += "Surface Loop(%d) = {%d,%d,%d,%d,%d,%d};\n"\
                %(self.ind_2d,self.ind_2d-1,self.ind_2d-3,self.ind_2d-5,self.ind_2d-7,self.ind_2d-9,self.ind_2d-11)
            self.ind_2d += 1
            self.file_string += "Volume(%d) = {%d};\n"\
                %(self.ind_2d,self.ind_2d-1)
            self.ind_vol = self.ind_2d
            self.ind_2d += 1
            self.file_string += "Physical Volume(%d) = {%d};\n"\
                %(self.ind_phys_vol, self.ind_2d-1)
        return

    def addSurface(self, p1, p2, p3, p4, scale, option):
        self.addPoint(p1, scale)
        self.addPoint(p2, scale)
        self.addPoint(p3, scale)
        self.addPoint(p4, scale)
        # print p1, p2, p3, p4
        l4 = self.addLineByIndex(self.ind_point-4, self.ind_point-3, scale)
        l3 = self.addLineByIndex(self.ind_point-3, self.ind_point-2, scale)
        l2 = self.addLineByIndex(self.ind_point-2, self.ind_point-1, scale)
        l1 = self.addLineByIndex(self.ind_point-1, self.ind_point-4, scale)

        self.addLineLoop(l4,l3,l2,l1,option)

    #creates the outer boundary, and sets it as a volume.
    def addOuterBound(self, p1, p2, scale, option):
        self.addBox(p1, p2, scale, option)
        # if option == "bound"
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
        # self.file_string += 'Background Field = %d;\n'%(self.ind_field)
        self.ind_field += 1
        return self.ind_field-1

    #Add a minimum field, based on the given list of fields. The minimum field applies the minimum resolution of all fields.
    def addMaxField(self, field_list):
        #create comma separated string of field indexes
        fields_string = ",".join(map(str, field_list))
        self.file_string += 'Field[%d] = Max;\n'%(self.ind_field)
        self.file_string += 'Field[%d].FieldsList = {%s};\n'%(self.ind_field, fields_string)
        # self.file_string += 'Background Field = %d;\n'%(self.ind_field)
        self.ind_field += 1
        return self.ind_field-1

    #Add a minimum field, based on the given list of fields. The minimum field applies the minimum resolution of all fields.
    def addMeanField(self, field_list, delta):
        #create comma separated string of field indexes
        fields_string = ",".join(map(str, field_list))
        self.file_string += 'Field[%d] = Mean;\n'%(self.ind_field)
        self.file_string += 'Field[%d].IField = %s;\n'%(self.ind_field, fields_string)
        self.file_string += 'Field[%d].Delta = %.15f;\n'%(self.ind_field, delta)
        # self.file_string += 'Background Field = %d;\n'%(self.ind_field)
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


# mw = MeshWriter()
# mw.addBox([0.0,0.0,0.0], [1.0,1.0,1.0], 1, "bound")
# mw.addBox([0.4,0.4,0.4], [0.6,0.6,0.6], 1, "seam")
# fields = []
# mw.addSurface([0.01,0.01,0.8],[0.99,0.01,0.8],[0.99,0.99,0.8],[0.01,0.99,0.8],1, "seam")
# fields += [mw.addTHField(0.5, 1, 0.01, 0.1)]
# mw.addSurface([0.0,0.0,0.5],[1.0,0.0,0.5],[1.0,1.0,0.5],[0.0,1.0,0.5],1, "bound")
# fields += [mw.addTHField(0.5, 1, 0.1, 0.3)]
# mw.addMinField(fields)
# with open('../data/domain.geo', 'w') as f: f.write(mw.file_string)
# subprocess.call(['gmsh -3 ../data/domain.geo'], shell=True)
# subprocess.call(['dolfin-convert ../data/domain.msh domain.xml'], shell=True)
# mesh = dolfin.Mesh('domain.xml')
# file_string = "../data/mesh.pvd"
# dolfin.File(file_string) << mesh
