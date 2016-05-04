#!/usr/bin/python

import networkx as nx
import matplotlib
import sys
import numpy as np

mcl_loops = int(sys.argv[2])

with open(sys.argv[1], 'r') as reader:
	G = nx.DiGraph()
	content = reader.readlines()

	for row in content:
		gene1, gene2, weight = row.split()
		G.add_edge(int(gene1), int(gene2), weight=abs(float(weight)))

#read through connected components
components = nx.weakly_connected_components(G)

#find clusters
for component in components:
	print component

	#create matrix
	adjacency_matrix = nx.to_numpy_matrix(G, nodelist=component)
	print adjacency_matrix
	cpy = adjacency_matrix

	#normalize matrix
	for j in xrange(0, len(adjacency_matrix)):
		sum = 0
		for k in xrange(len(adjacency_matrix)):
			sum += adjacency_matrix.item((k,j))
		for m in xrange(len(adjacency_matrix)):
			if sum != 0:
				cpy[m,j] = adjacency_matrix.item(m,j)/sum
	#print cpy
	old = cpy

	for l in xrange(mcl_loops):
		exp = np.dot(old, old)
		#print exp

		#inflation:
		#square all elements in matrix
		inf = np.square(exp)
		#print inf

		#normalize matrix
		new = inf
		for j in range(0, len(inf)):
			sum = 0
			for k in range(0, len(inf)):
				sum = sum + inf.item((k,j))
			for m in range(0, len(inf)):
				if sum != 0:
					new[m,j] = inf.item(m,j)/sum
					if new[m,j] < 0.000001:
						new[m,j] = 0
		old = new
	print "Converged to: \n", old
	#find cluster!
	#every pair (i, j) which has a non-zero value are in the same cluster

	#calculate change of matrices
	#diff = np.subtract(new, old)
	#print diff
