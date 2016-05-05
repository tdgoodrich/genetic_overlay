#!/usr/bin/python

import networkx as nx
import matplotlib
import sys
import csv
import numpy as np
import time


def column_stochastic_matrix(matrix):
	# Row-stochasticize matrix
	for col in xrange(len(matrix)):
		col_sum = np.sum(matrix[:,col])
		if col_sum != 0:
			matrix[:,col] = np.divide(matrix[:,col], col_sum)
	return matrix

def threshold_reduction(element, threshold=0.00001):
	if element < threshold:
		return 0
	else:
		return element

threshold_reduction_vectorized = np.vectorize(threshold_reduction)

# Note: r=1.4 chosen to match original study
def markov_clustering(filename, max_loops=20, r=1.4):
	# Read adjacency matrix
	with open(sys.argv[1], 'r') as infile:
		data = list(csv.reader(infile))
		genes = [x.replace(" ", "") for x in data[0]]

		# Remove headers
		del data[0]
		del genes[0]

		for i in xrange(len(data)):
			row = data[i]
			del row[0] # Remove header
			data[i] = [abs(float(v)) for v in row]

	# Construct the graph
	adjacency_matrix = np.array(data)
	G = nx.from_numpy_matrix(adjacency_matrix, create_using=nx.DiGraph())

	# Loop over each component
	components = nx.weakly_connected_components(G)
	for component in components:
		# Create component's matrix
		component_matrix = nx.to_numpy_matrix(G, nodelist=component)
		#print "Component's adjacency matrix: \n", component_matrix

		# Stochastcize the columns
		component_matrix = column_stochastic_matrix(component_matrix)

		for l in xrange(max_loops):
			# Expansion step
			component_matrix = np.dot(component_matrix, component_matrix)

			# Inflation step
			# Assumes we want to square
			component_matrix = np.power(component_matrix, r)

			# Normalize matrix and threshold close-to-zero
			component_matrix = column_stochastic_matrix(component_matrix)
			component_matrix = threshold_reduction_vectorized(component_matrix)

		#print component_matrix

		# Assumes n by n matrix
		clusters = set()
		n = len(component_matrix)
		for row in xrange(n):
			attracted = [i for i in xrange(n) if component_matrix[row,i] > 0]
			if attracted != []:
				cluster = frozenset(attracted + [row])
				clusters.add(cluster)

		# Build lookup table for a gene's clusters
		#print "Genes: ", genes
		cluster_lookup = {}
		cluster_id = 1
		for cluster in clusters:
			#print cluster
			for element in cluster:
				cluster_lookup[genes[element]] = cluster_lookup.get(genes[element], []) + [cluster_id]
			cluster_id += 1
		return cluster_lookup

def write_cluster_lookup(cluster_lookup, filename):
	with open(filename, "w") as outfile:
		for gene in cluster_lookup:
			for cluster in cluster_lookup[gene]:
				outfile.write("%s %s\n" % (gene, cluster))

start = time.time()
cluster_lookup = markov_clustering(sys.argv[1], int(sys.argv[2]))
stop = time.time()
write_cluster_lookup(cluster_lookup, "Data/mcl_clusters.txt")
print "Total run time: %.4f sec" % (stop - start)
