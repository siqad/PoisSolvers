
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

def point_inside_polygon(x, y, poly, include_edges=True):
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
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]
    for i in range(1, n + 1):
        p2x, p2y = poly[i % n]
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
