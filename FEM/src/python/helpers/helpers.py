 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Convenience functions for general data manipulation

import numpy as np
import itertools
#elec_list and elec_poly_list are a list of electrodes and list of polygonal electrodes
def getBB(sqconn):
    elec_list = getElectrodeCollections(sqconn)
    sim_params = sqconn.getAllParameters()
    bcs = sim_params["bcs"]
    padding = float(sim_params["padding"])
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
    if padding == 0.0:
        if bcs == "periodic":
            scale = 0.2
        else:
            scale = 4.0
        xs = [min_x-scale*(max_x-min_x), max_x+scale*(max_x-min_x)]
        ys = [min_y-scale*(max_y-min_y), max_y+scale*(max_y-min_y)]
        return xs, ys
    else:
        xs = [min_x-abs(padding), max_x+abs(padding)]
        ys = [min_y-abs(padding), max_y+abs(padding)]
        return xs, ys

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
def adjustBoundaries(xs, ys, m_p):
    candidates = np.array([])
    for key in m_p:
        if m_p[key]:
            candidates = np.append(candidates, m_p[key][0])
            candidates = np.append(candidates, sum(m_p[key]))
    zmin = -np.max(np.abs(candidates))*1.5
    zmax = -zmin
    b_di = 0.0 #at the surface.
    return xs[0], xs[1], ys[0], ys[1], zmin, zmax, b_di
