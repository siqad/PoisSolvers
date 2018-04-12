import pygmsh as pg

geom = pg.built_in.Geometry()
shape = geom.add_box(0, 1, 0, 1, 0, 1, 0.05)

# print geom
points, cells, point_data, cell_data, field_data = pg.generate_mesh(geom)

# print points
# print cells