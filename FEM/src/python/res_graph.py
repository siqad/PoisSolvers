import networkx as nx
import itertools
import matplotlib.pyplot as plt
import shapely.geometry as sp
import numpy as np

class ResGraph():

    def __init__(self, elec_dict, elec_list, dir):
        self.g = nx.Graph()
        self.elec_dict = elec_dict
        self.elec_list = elec_list
        self.dir = dir
        self.z_bounds = []
        self.edge_ind = 0
        self.node_ind = 0
        self.setZMaxes()

        self.buildGraph()

        self.debugPrint()

    def setZMaxes(self):
        self.z_bounds.append(min(elec.z1 for elec in self.elec_list))
        self.z_bounds.append(max(elec.z2 for elec in self.elec_list))
        # print(self.z_bounds)

    def addFromOverlap(self, pair):
        a, b = pair

        #if any of these are true, no overlap
        if a.x1 > b.x2 or a.x2 < b.x1 \
            or a.y1 > b.y2 or a.y2 < b.y1 \
            or a.z1 > b.z2 or a.z2 < b.z1:
            pass

        #Overlap exists between electrodes, this pair is connected.
        #Find which coordinate overlaps, ignore rotated ones for now.
        else:
            #y and z coordinates form the intersection
            if a.x1 == b.x2:
                box_a = sp.box(a.y1, a.z1, a.y2, a.z2)
                box_b = sp.box(b.y1, b.z1, b.y2, b.z2)
                inter = box_a.intersection(box_b)
                pos = [a.x1, inter.centroid.x, inter.centroid.y]
            if a.x2 == b.x1:
                box_a = sp.box(a.y1, a.z1, a.y2, a.z2)
                box_b = sp.box(b.y1, b.z1, b.y2, b.z2)
                inter = box_a.intersection(box_b)
                pos = [a.x2, inter.centroid.x, inter.centroid.y]

            # x and z coordinates form the intersection
            if a.y1 == b.y2:
                box_a = sp.box(a.x1, a.z1, a.x2, a.z2)
                box_b = sp.box(b.x1, b.z1, b.x2, b.z2)
                inter = box_a.intersection(box_b)
                pos = [inter.centroid.x, a.y1, inter.centroid.y]
            if a.y2 == b.y1:
                box_a = sp.box(a.x1, a.z1, a.x2, a.z2)
                box_b = sp.box(b.x1, b.z1, b.x2, b.z2)
                inter = box_a.intersection(box_b)
                pos = [inter.centroid.x, a.y2, inter.centroid.y]

            #x and y coordinates form the intersection
            if a.z1 == b.z2:
                box_a = sp.box(a.x1, a.y1, a.x2, a.y2)
                box_b = sp.box(b.x1, b.y1, b.x2, b.y2)
                inter = box_a.intersection(box_b)
                pos = [inter.centroid.x, inter.centroid.x, a.z1]
            if a.z2 == b.z1:
                box_a = sp.box(a.x1, a.y1, a.x2, a.y2)
                box_b = sp.box(b.x1, b.y1, b.x2, b.y2)
                inter = box_a.intersection(box_b)
                pos = [inter.centroid.x, inter.centroid.x, a.z2]
            includes = [a.id, b.id]
            curr_node = self.addNode(pos=pos, includes=includes, end=False)

    # Add in the nodes that connect to the ceiling
    def addCeilingNodes(self, net_elec):
        for item in net_elec:
            #check if electrode connects to the ceiling
            if item.z2 == self.z_bounds[1]:
                #find the center of the intersection (in x and y)
                box = sp.box(item.x1, item.y1, item.x2, item.y2)
                center = box.centroid
                #add its center
                curr_node = self.addNode(pos=[center.x, center.y,item.z2], includes=[item.id], end=False)
                self.addEdge(curr_node, "ceiling")

    def addFloorNodes(self, net_elec):
        for item in net_elec:
            #check if electrode connects to the floor
            if item.z1 == self.z_bounds[0]:
                #find the center of the intersection (in x and y)
                box = sp.box(item.x1, item.y1, item.x2, item.y2)
                center = box.centroid
                #add its center
                curr_node = self.addNode(pos=[center.x, center.y,item.z1], includes=[item.id], end=False)
                self.addEdge(curr_node, "floor")

    # Adds an edge to the graph between two nodes a and b, where a and b are node keys.
    def addEdge(self, a, b, **kwargs):
        self.g.add_edge(a, b, **kwargs)

    # Adds a node to the graph at with a specified key. If no key is specified, the internal node index is used and incremented.
    def addNode(self, key=None, **kwargs):
        if key == None:
            key = self.node_ind
            self.node_ind += 1
        self.g.add_node(key, **kwargs)
        return key

    #nodes defined by electrode intersections
    def buildNodes(self, key):
        # Create all possible pair combinations of electrodes within a net
        pairs = itertools.combinations(self.elec_dict[key], 2)
        # Check for overlap
        for pair in pairs:
            self.addFromOverlap(pair)

    #gets the nodes that the electrode with a given id is involved with
    def getRelevantNodes(self, id):
        relevant_nodes = []
        for node in self.g.nodes(data='includes'):
            if node[1] != None and id in node[1]:
                relevant_nodes.append(node[0])
        return relevant_nodes

    def getDistanceDict(self, subgraph, rel_nodes):
        node_pos = subgraph.nodes(data='pos')
        pairs = [x for x in itertools.combinations(rel_nodes, 2) ]
        distances = []

        for pair in pairs:
            a, b = pair
            point_a = np.array(node_pos[a])
            point_b = np.array(node_pos[b])
            #Get the distance between node positions
            distances.append(np.linalg.norm(point_a - point_b))
        return dict(zip(pairs, distances))

    # Connect each node to its two nearest neighbours. End nodes are connected only to a single nearest neighbour.
    def connectSubgraph(self, sg, dist_dict):
        #once for each node
        for node, attributes in sg.nodes(data=True):
            #sort the pairs based on their distances
            sorted_pairs = sorted(dist_dict, key=lambda x: (dist_dict[x]))
            #get the pairs relevant to this node
            sorted_filtered_pairs = list(filter(lambda x: node in x, sorted_pairs))

            for pair in sorted_filtered_pairs:
                #if a path exists between these two nodes, try next pair
                if nx.has_path(sg, pair[0],pair[1]):
                    continue
                #no path exists between this pair of closest nodes. Connect them!
                a, b = pair
                self.addEdge(a,b)
                print("Connecting ", a, b)
                #this node is now connected in the graph. head to the next node.
                break

            # print(sorted_filtered_pairs)

            #for each entry in the dist
            # for pair in dist_dict:
            #     #Only about the pair if it includes the current node of interest
            #     if node in pair:
            #         print(pair, node)
            #         #get the relevant pair keys in ascending order.
            #         sorted_pairs = sorted(dist_dict, key=lambda x: (dist_dict[x]))
            #         print(sorted_pairs)
            #         sorted_filtered_pairs = list(filter(lambda x: node in x, sorted_pairs))
            #         if attributes["end"] == True:
            #             #End item, connect only the closest one.
            #             a, b = sorted_filtered_pairs[0]
            #             self.addEdge(a, b)
            #             print("Connecting: ", a, b)
            #         else:
            #             #Connect the two closest
            #             for i in range(2):
            #                 a, b = sorted_filtered_pairs[i]
            #                 self.addEdge(a, b)
            #                 print("Connecting: ", a, b)

    #mark the two nodes farthest from each other as being the ends
    def markMaxDist(self, dist_dict, sg):
        max_key = max(dist_dict, key=lambda x: dist_dict[x])
        # max_key should be a list with two node keys
        for node in sg.nodes():
            # print("Node: ", node, "max_key: ", max_key)
            if node in max_key:
                self.addNode(node, end=True)
                # print("Setting node ", node, " as end node.")

    #edges are defined by the electrode geometries
    def buildEdges(self):
        for item in self.elec_list:
            rel_nodes = self.getRelevantNodes(item.id)
            sg = self.g.subgraph(rel_nodes)
            dist_dict = self.getDistanceDict(sg, rel_nodes)
            # print(dist_dict)
            self.markMaxDist(dist_dict, sg)
            self.connectSubgraph(sg, dist_dict)
            #reset the end attribute for all nodes
            nx.set_node_attributes(sg, False, "end")

    def buildGraph(self):
        self.g = nx.Graph()

        self.addNode("ceiling")
        self.addNode("floor")

        #Once per net
        for key in self.elec_dict.dict:
            #identify the electrodes that have the highest height.
            self.addCeilingNodes(self.elec_dict[key])
            #identify the electrodes that have the lowest height.
            self.addFloorNodes(self.elec_dict[key])
            #build intermediate nodes
            self.buildNodes(key)
            self.buildEdges()

        self.exportGraph()

    def exportGraph(self):
        plt.figure()
        nx.draw(self.g, with_labels=True)
        plt.savefig(self.dir+"/graph.pdf", bbox_inches='tight')

    def debugPrint(self):
        print(list(self.g.nodes(data=True)))

# if __name__ == "__main__":
#     pass
