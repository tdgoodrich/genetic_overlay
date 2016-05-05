#!/usr/bin/python

import sys
import networkx as nx
import data_processing as dp

################################################################################
# Main
################################################################################

if __name__ == "__main__":
    if sys.argv[1] == "build_clean_data":
        vertices, edgelist = dp.build_threshold_data(threshold=0.08)
        dp.write_edgelist(edgelist)
        dp.write_vertices(vertices)
        dp.write_adjacency_matrices(vertices, edgelist)
    elif sys.argv[1] == "check_arcs":
        dp.check_arcs()
    elif sys.argv[1] == "build_edgelist":
        vertices, edgelist = dp.build_edgelist()
        dp.write_edgelist(edgelist)
        dp.write_vertices(vertices)
        dp.write_adjacency_matrices(vertices, edgelist)
    elif sys.argv[1] == "plot_graph":
        vertices, edgelist = dp.build_graph()
        dp.plot_graph(vertices, edgelist)
    elif sys.argv[1] == "gephi":
        G = dp.build_graph()
        nx.write_gexf(G, "Data/network.gexf")
