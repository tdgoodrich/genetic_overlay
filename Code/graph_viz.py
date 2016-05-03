#!/usr/bin/python

import networkx as nx
import matplotlib
import sys
import matplotlib.pyplot as plt


with open(sys.argv[1], 'r') as reader:

	G=nx.Graph()

	content = reader.readlines()

	for row in content:
		vals = row.split()
		gene1 = vals[0]
		gene2 = vals[1]
		weight = vals[2]
		G.add_edge(gene1, gene2)

	print G.number_of_nodes()

	print G.number_of_edges()

	pos = nx.spring_layout(G)

	nx.draw(G, pos)

	plt.savefig("graph.png")