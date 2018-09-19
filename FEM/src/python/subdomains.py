 # @author:   Nathan
 # @created:  2018.08.24
 # @editted:  2017.08.24 - Nathan
 # @license:  GNU LGPL v3
 #
 # @desc:     Reimplementations for the boundaries and subdomains

import dolfin

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




# this is correct
# class PeriodicBoundary(dolfin.SubDomain):
#     def __init__(self, x_min, x_max, y_min, y_max):
#         self.x_min = x_min
#         self.x_max = x_max
#         self.y_min = y_min
#         self.y_max = y_max
#         dolfin.SubDomain.__init__(self)
#
#     def inside(self, x, on_boundary):
#         # the boundary minus the top and right sides
#         return on_boundary and not (dolfin.near(x[0], self.x_max) or dolfin.near(x[1], self.y_max))
#     def map(self, x, y):
#             if dolfin.near(x[0], self.x_max):
#                 y[0] = self.x_min
#             elif dolfin.near(x[1], self.y_max):
#                 y[1] = self.x_min
#             else:
#                 y[i] = x[i]


# Sub domain for Periodic boundary condition in x and y
class PeriodicBoundary(dolfin.SubDomain):
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        print(x_min, x_max, y_min, y_max)
        dolfin.SubDomain.__init__(self)
    def inside(self, x, on_boundary):
        return bool((dolfin.near(x[0], self.x_min) or dolfin.near(x[1], self.y_min)) and
                (not ((dolfin.near(x[0], self.x_min) and dolfin.near(x[1], self.y_max)) or
                        (dolfin.near(x[0], self.x_max) and dolfin.near(x[1], self.y_min)))) and on_boundary)
    def map(self, x, y):
        if dolfin.near(x[0], self.x_max) and dolfin.near(x[1], self.y_max):
            y[0] = self.x_min
            y[1] = self.y_min
            y[2] = x[2]
        elif dolfin.near(x[0], self.x_max):
            y[0] = self.x_min
            y[1] = x[1]
            y[2] = x[2]
        elif dolfin.near(x[1], self.y_max):
            y[0] = x[0]
            y[1] = self.y_min
            y[2] = x[2]
        else:
            y[0] = 1000
            y[1] = 1000
            y[2] = 1000
    # def map(self, x, y):
    #     if dolfin.near(x[0], self.x_max) and dolfin.near(x[1], self.y_max):
    #         y[0] = x[0] - self.x_max
    #         y[1] = x[1] - self.y_max
    #         y[2] = x[2]
    #     elif dolfin.near(x[0], self.x_max):
    #         y[0] = x[0] - self.x_max
    #         y[1] = x[1]
    #         y[2] = x[2]
    #     elif dolfin.near(x[1], self.y_max):
    #         y[0] = x[0]
    #         y[1] = x[1] - self.y_max
    #         y[2] = x[2]
    #     else:
    #         y[0] = 1000
    #         y[1] = 1000
    #         y[2] = 1000
        # if dolfin.dolfin.near(x[0], self.x_max) and dolfin.dolfin.near(x[1], self.y_max):
        #     y[0] = x[0] - self.x_max
        #     y[1] = x[1]
        #     y[2] = x[2]
        #     # y[0] = x[0] - self.x_max
        #     # y[1] = x[1] - self.y_max
        #     # y[2] = x[2]
        # elif dolfin.near(x[0], self.x_max):
        #     y[0] = x[0] - self.x_max
        #     y[1] = x[1]
        #     y[2] = x[2]
        # else:   # near(x[1], 1)
        #     y[0] = x[0]
        #     y[1] = x[1] - self.y_max
        #     y[2] = x[2]

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
    def __init__(self, xs, ys, zs):
        self.xs = xs
        self.ys = ys
        self.zs = zs
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return (dolfin.between(x[0], (self.xs[0], self.xs[1])) \
            and dolfin.between(x[1], (self.ys[0], self.ys[1])) \
            and dolfin.between(x[2], (self.zs[0], self.zs[1])) )

# INTERNAL BOUNDARY CONDITION ELECTRODEPOLY
class ElectrodePoly(dolfin.SubDomain):
    def __init__(self, vertices, zs):
        self.vertices = vertices
        self.zs = zs
        dolfin.SubDomain.__init__(self) # Call base class constructor!
    def inside(self, x, on_boundary):
        return ( self.point_inside_polygon(x[0], x[1]) \
            and dolfin.between(x[2], (self.zs[0], self.zs[1])) )
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
