import os
import numpy as np
import json

class Exporter():
    def __init__(self, in_path="", out_path="", json_path=None, sqconn=None):
        print("Exporter created")
        if in_path.strip():
            self.in_path = in_path
            self.abs_in_dir = os.path.abspath(os.path.dirname(in_path))
        if out_path.strip():
            self.out_path = out_path
            self.abs_out_dir = os.path.abspath(os.path.dirname(out_path))

        self.json_export_path = json_path
        self.sqconn = sqconn
        self.db_hist = []
        self.db_pots_export = {}
        self.db_pots_export['pots'] = []

    def exportDBs(self, step, steps, db_list, u, dielec):
        if db_list:
            db_pots = []
            pots_snapshot = []
            for db in db_list:
                db_pots.append([2*np.pi*step/steps, db.x, db.y, u(db.x, db.y, dielec)])
                pots_snapshot.append(u(db.x, db.y, dielec))
            self.db_hist.extend(db_pots)
            self.db_pots_export['pots'].append(pots_snapshot)

    def exportPotentialXML(self, data_2d):
        X, Y, Z, nx, ny = data_2d
        XYZ = []
        for i in range(nx):
            for j in range(ny):
                XYZ.append([X[i,j],Y[i,j],Z[i,j]])
        self.sqconn.export(potential=XYZ)

    def exportDBHistory(self):
        self.sqconn.export(db_pot=self.db_hist)

    def exportDBJSON(self):
        if self.json_export_path and len(self.db_pots_export['pots']) > 0:
            with open(self.json_export_path, 'w') as export_file:
                json.dump(self.db_pots_export, export_file)
