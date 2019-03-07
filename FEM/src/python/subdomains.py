 # @author:   Nathan
 # @created:  2018.08.24
 # @editted:  2017.08.24 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Reimplementations for the boundaries and subdomains

import dolfin
import numpy as np
import itertools
# Create classes for defining parts of the boundaries and the interior
# of the domain
class Left(dolfin.SubDomain): #x_min
    def __init__(self, boundary_x_min):
        self.boundary_x_min = boundary_x_min
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], self.boundary_x_min)

class Right(dolfin.SubDomain): #x_max
    def __init__(self, boundary_x_max):
        self.boundary_x_max = boundary_x_max
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[0], self.boundary_x_max)

class Bottom(dolfin.SubDomain): #y_min
    def __init__(self, boundary_y_min):
        self.boundary_y_min = boundary_y_min
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], self.boundary_y_min)

class Top(dolfin.SubDomain): #y_max
    def __init__(self, boundary_y_max):
        self.boundary_y_max = boundary_y_max
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[1], self.boundary_y_max)

class Back(dolfin.SubDomain): #z_min
    def __init__(self, boundary_z_min):
        self.boundary_z_min = boundary_z_min
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], self.boundary_z_min)

class Front(dolfin.SubDomain): #z_max
    def __init__(self, boundary_z_max):
        self.boundary_z_max = boundary_z_max
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return dolfin.near(x[2], self.boundary_z_max)

# Sub domain for Periodic boundary condition in x and y
class PeriodicBoundary(dolfin.SubDomain):
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = float(x_min)
        self.x_max = float(x_max)
        self.y_min = float(y_min)
        self.y_max = float(y_max)
        # print("Periodic:", self.x_min, self.x_max, self.y_min, self.y_max)
        # print(self.inside([-200,-100,10], True))
        dolfin.SubDomain.__init__(self, map_tol=1e-2)
        print(self.map_tolerance)
    def inside(self, x, on_boundary):
        return bool((dolfin.near(x[0], self.x_min) or dolfin.near(x[1], self.y_min)) and
                (not ((dolfin.near(x[0], self.x_min) and dolfin.near(x[1], self.y_max)) or
                        (dolfin.near(x[0], self.x_max) and dolfin.near(x[1], self.y_min)))) and on_boundary)

    def map(self, x, y):
        if dolfin.near(x[0], self.x_max) and dolfin.near(x[1], self.y_max):
            y[0] = x[0] - self.x_max + self.x_min
            y[1] = x[1] - self.y_max + self.y_min
            # y[0] = self.x_min
            # y[1] = self.y_min
            y[2] = x[2]
        # elif dolfin.near(x[0], self.x_max, 1E-3):
        elif dolfin.near(x[0], self.x_max):
            # print("x: ", x[0] - (self.x_max - self.x_min), x[1], x[2])
            # y[0] = x[0]
            # y[0] = x[0] - (self.x_max - self.x_min)
            # print("x: ", x[0], x[0] - self.x_max + self.x_min)
            # print("x: ", x[0], x[0] - self.x_max + self.x_min)
            y[0] = x[0] - self.x_max + self.x_min
            y[1] = x[1]
            y[2] = x[2]
            # print(self.x_min, y[0], x[0])
        elif dolfin.near(x[1], self.y_max):
            # print("y: ", x[1], x[1] - self.y_max + self.y_min)
            # print("y: ", x[1] - (self.y_max - self.y_min))
            y[0] = x[0]
            # y[1] = x[1] - (self.y_max - self.y_min)
            y[1] = x[1] - self.y_max + self.y_min
            # y[1] = self.y_min
            y[2] = x[2]
        else:
            y[0] = (self.x_max + self.x_min)/2
            y[1] = (self.y_max + self.y_min)/2
            y[2] = x[2]

# INTERNAL BOUNDARY CONDITION FOR DIELECTRIC
class Air(dolfin.SubDomain):
    def __init__(self, bounds):
        self.boundary_dielectric = bounds[0]
        self.boundary_z_max = bounds[1]
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return (dolfin.between(x[2], (self.boundary_dielectric, self.boundary_z_max)))

# INTERNAL BOUNDARY CONDITION ELECTRODE
class Electrode(dolfin.SubDomain):
    def __init__(self, electrode):
        # extract the dimensions from the electrode
        self.electrode = electrode
        self.tol = 1e-2
        self.xs = [electrode.x1, electrode.x2]
        self.ys = [electrode.y1, electrode.y2]
        self.zs = [electrode.z1, electrode.z2]
        self.vertices = []
        self.getVertices()
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def getVertices(self):
        #translate the vertices to the origin
        theta = np.deg2rad(self.electrode.angle)
        origin = np.array([sum(self.xs)/2, sum(self.ys)/2])
        centered_x = [x - origin[0] for x in self.xs]
        centered_y = [y - origin[1] for y in self.ys]
        for vertex in itertools.product(centered_x, centered_y):
        #     #for each vertex rotate about the origin
            new_vertex = np.array([vertex[0]*np.cos(theta)-np.sin(theta)*vertex[1], \
                          vertex[0]*np.sin(theta)+vertex[1]*np.cos(theta)])
            self.vertices.append(tuple(new_vertex + origin))
        # print(self.vertices)
            # print(tuple(new_vertex + origin))
    def inside(self, x, on_boundary):
        if np.abs(self.electrode.angle) > 1e-16:
            if self.point_inside_polygon(x[0],x[1]):
                # print(x, "True without modification.")
                return True
            elif self.point_inside_polygon(x[0]+self.tol,x[1]) or \
                 self.point_inside_polygon(x[0]+self.tol,x[1]+self.tol) or \
                 self.point_inside_polygon(x[0]+self.tol,x[1]-self.tol) or \
                 self.point_inside_polygon(x[0]-self.tol,x[1]) or \
                 self.point_inside_polygon(x[0]-self.tol,x[1]+self.tol) or \
                 self.point_inside_polygon(x[0]-self.tol,x[1]-self.tol) or \
                 self.point_inside_polygon(x[0],x[1]+self.tol) or \
                 self.point_inside_polygon(x[0],x[1]-self.tol):
                # print(x, "True after modification.")
                return True
            else:
                # print(x, "Outside")
                return False
        else:
            return ( self.point_inside_polygon(x[0], x[1]) \
                and dolfin.between(x[2], (self.zs[0], self.zs[1])))

    def point_inside_polygon(self, x, y, include_edges=True):
        '''
        Test if point (x,y) is inside polygon poly.

        poly is N-vertices polygon defined as
        [(x1,y1),...,(xN,yN)] or [(x1,y1),...,(xN,yN),(x1,y1)]
        (function works fine in both cases)

        Geometrical idea: point is inside polygon if horisontal beam
        to the right from point crosses polygon even number of times.
        Works fine for non-convex polygons.

        Returns True if point (x, y) is inside poly, else False
        '''
        n = len(self.vertices)
        inside = False
        p1x, p1y = self.vertices[0]
        for i in range(1, n + 1):
            p2x, p2y = self.vertices[i % n]
            if p1y == p2y:
                if y == p1y:
                    if min(p1x, p2x) <= x <= max(p1x, p2x):
                        # point is on horizontal edge
                        inside = include_edges
                        break
                    elif x < min(p1x, p2x):  # point is to the left from current edge
                        inside = not inside
            else:  # p1y!= p2y
                if min(p1y, p2y) <= y <= max(p1y, p2y):
                    xinters = (y - p1y) * (p2x - p1x) / float(p2y - p1y) + p1x
                    if x == xinters:  # point is right on the edge
                        inside = include_edges
                        break
                    if x < xinters:  # point is to the left from current edge
                        inside = not inside
            p1x, p1y = p2x, p2y
        return inside
