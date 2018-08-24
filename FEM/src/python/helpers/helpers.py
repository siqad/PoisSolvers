 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Convenience functions for general data manipulation

import numpy as np

#elec_list and elec_poly_list are a list of electrodes and list of polygonal electrodes
def getBB(elec_list, elec_poly_list):
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
    xs = [min_x-2.0*(max_x-min_x), max_x+2.0*(max_x-min_x)]
    ys = [min_y-2.0*(max_y-min_y), max_y+2.0*(max_y-min_y)]
    return xs, ys

def getMetalParams(sqconn):
    for layer in sqconn.getLayers():
        if layer.name == "Metal":
            metal_offset = float(layer.zoffset)
            metal_thickness = float(layer.zheight)
            break
    return metal_offset, metal_thickness

def getElectrodeCollections(sqconn):
    m_per_A = 1.0E-10 #metres per angstrom
    elec_list = []
    for elec in sqconn.electrodeCollection():
        elec_curr = elec
        #units given are in angstroms, convert to metres
        elec_curr.x1 *= m_per_A
        elec_curr.x2 *= m_per_A
        elec_curr.y1 *= m_per_A
        elec_curr.y2 *= m_per_A
        elec_list.append(elec_curr)
    return elec_list

def getElectrodePolyCollections(sqconn):
    m_per_A = 1.0E-10 #metres per angstrom
    elec_poly_list = []
    for elec_poly in sqconn.electrodePolyCollection():
        # pass
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

def adjustBoundaries(xmin, xmax, ymin, ymax, moff, mthick):
    if xmin == 0:
        xmin -= 0.01*xmax
    if ymin == 0:
        ymin -= 0.01*ymax
    zmin = -np.max(np.array([np.abs(moff), np.abs(mthick)]))*10.0
    zmax = -zmin
    b_di = 0.0 #at the surface.
    return xmin, xmax, ymin, ymax, zmin, zmax, b_di
