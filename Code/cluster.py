#!/usr/bin/python

import networkx as nx
import matplotlib
import sys
import numpy as np

mcl_loops = int(sys.argv[2])

with open(sys.argv[1], 'r') as reader:

	G=nx.Graph()

	content = reader.readlines()

	for row in content:
		vals = row.split()
		gene1 = vals[0]
		gene2 = vals[1]
		wt = float(vals[2])
		wt = abs(wt)
		G.add_edge(gene1, gene2, weight=wt)

#A = nx.adjacency_matrix(G)
#print(A.todense())

#read through connected components
ccs = nx.connected_components(G)

#find clusters
for i in ccs:
	#create matrix
	ccm = nx.to_numpy_matrix(G, nodelist=i)
	cpy = ccm
	#normalize matrix
	for j in range(0, len(ccm)):
		sum = 0
		for k in range(0, len(ccm)):
			sum = sum + ccm.item((k,j))
		for m in range(0, len(ccm)):
			if sum != 0:
				cpy[m,j] = ccm.item(m,j)/sum
	#print cpy
	old = cpy

	for l in xrange(0,mcl_loops):
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
				sum = sum + ccm.item((k,j))
			for m in range(0, len(inf)):
				if sum != 0:
					new[m,j] = inf.item(m,j)/sum
		old = new
	#print old
	#find cluster!
	#every pair (i, j) which has a non-zero value are in the same cluster

	#calculate change of matrices
	#diff = np.subtract(new, old)
	#print diff