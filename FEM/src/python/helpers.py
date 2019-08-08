 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Convenience functions for general data manipulation

import numpy as np
import itertools
#elec_list and elec_poly_list are a list of electrodes and list of polygonal electrodes
def getBB(sqconn, pad):
    elec_list = getElectrodeCollections(sqconn)
    sim_params = sqconn.getAllParameters()
    x_list = []
    y_list = []

    if elec_list:
        for elec in elec_list:
            #translate the verticecs to the origin
            theta = np.deg2rad(elec.angle)
            xs = [elec.x1,elec.x2]
            ys = [elec.y1,elec.y2]
            vertices = []
            origin = np.array([sum(xs)/2, sum(ys)/2])
            centered_x = [x - origin[0] for x in xs]
            centered_y = [y - origin[1] for y in ys]
            for vertex in itertools.product(centered_x, centered_y):
            #     #for each vertex rotate about the origin
                new_vertex = np.array([vertex[0]*np.cos(theta)-np.sin(theta)*vertex[1], \
                              vertex[0]*np.sin(theta)+vertex[1]*np.cos(theta)])
                vertices.append(tuple(new_vertex + origin))
            x_list += [vertex[0] for vertex in vertices]
            y_list += [vertex[1] for vertex in vertices]

    min_x = min(x_list)
    max_x = max(x_list)
    min_y = min(y_list)
    max_y = max(y_list)
    xs = [min_x-abs(pad[0]), max_x+abs(pad[0])]
    ys = [min_y-abs(pad[1]), max_y+abs(pad[1])]
    return xs, ys, [min_x, max_x], [min_y, max_y]

def getMetalParams(sqconn):
    metal_params = {}
    layer_id = 0
    for layer in sqconn.getLayers():
        if layer.type == "Electrode":
            # Save layer ID and offsets
            metal_params[layer_id] = (float(layer.zoffset), float(layer.zheight))
        layer_id += 1
    return metal_params

def getElectrodeCollections(sqconn):
    elec_list = []
    for elec in sqconn.electrodeCollection():
        elec_list.append(elec)
    return elec_list

def getDBCollections(sqconn):
    db_list = []
    for db in sqconn.dbCollection():
        db_list.append(db)
    return db_list

# def adjustBoundaries(xmin, xmax, ymin, ymax, m_p):
def adjustBoundaries(xs, ys, m_p, gp_depth):
    candidates = np.array([])
    for key in m_p:
        if m_p[key]:
            candidates = np.append(candidates, m_p[key][0])
            candidates = np.append(candidates, sum(m_p[key]))
    zmax = np.max(np.abs(candidates))*1.5
    zmin = -gp_depth
    b_di = 0.0 #at the surface.
    return xs[0], xs[1], ys[0], ys[1], zmin, zmax, b_di
