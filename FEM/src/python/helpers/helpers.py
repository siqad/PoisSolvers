 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Convenience functions for general data manipulation

import numpy as np

#elec_list and elec_poly_list are a list of electrodes and list of polygonal electrodes
def getBB(elec_list, elec_poly_list, bcs, padding):
    x_list = []
    y_list = []
    if elec_list:
        x_list += [a.x1 for a in elec_list] + [b.x2 for b in elec_list]
        y_list += [a.y1 for a in elec_list] + [b.y2 for b in elec_list]
    if elec_poly_list:
        for elec_poly in elec_poly_list:
            x_list += [c[0] for c in elec_poly.vertex_list]
            y_list += [c[1] for c in elec_poly.vertex_list]
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
        if "Metal" in layer.name:
            # Save layer ID and offsets
            metal_params[layer_id] = (float(layer.zoffset), float(layer.zheight))
        layer_id += 1
    return metal_params

def getElectrodeCollections(sqconn):
    elec_list = []
    for elec in sqconn.electrodeCollection():
        elec_curr = elec
        elec_list.append(elec_curr)
    return elec_list

def getElectrodePolyCollections(sqconn):
    m_per_A = 1.0E-10 #metres per angstrom
    elec_poly_list = []
    for elec_poly in sqconn.electrodePolyCollection():
        elec_poly_curr = elec_poly
        #convert units to metres
        vertex_list = [list(n) for n in elec_poly_curr.vertices]
        for i in range(len(vertex_list)):
            vertex_list[i][0] *= m_per_A
            vertex_list[i][1] *= m_per_A
        elec_poly_curr.vertex_list = vertex_list
        elec_poly_list.append(elec_poly_curr)
    return elec_poly_list

def getDBCollections(sqconn):
    db_list = []
    for db in sqconn.dbCollection():
        db_list.append(db)
    return db_list

def adjustBoundaries(xmin, xmax, ymin, ymax, m_p):
    candidates = np.array([])
    for key in m_p:
        if m_p[key]:
            candidates = np.append(candidates, m_p[key][0])
            candidates = np.append(candidates, sum(m_p[key]))
    zmin = -np.max(np.abs(candidates))*1.5
    zmax = -zmin
    b_di = 0.0 #at the surface.
    return xmin, xmax, ymin, ymax, zmin, zmax, b_di
