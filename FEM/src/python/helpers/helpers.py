 # @author:   Nathan
 # @created:  2018.08.23
 # @editted:  2018.08.23 - Nathan
 # @license:  Apache License 2.0
 #
 # @desc:     Functions that create the GMSH .geo file

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
