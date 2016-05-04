#!/usr/bin/python

import networkx as nx
import matplotlib
import sys
import csv
import numpy as np

mcl_loops = int(sys.argv[2])
clusterList = []

#with open(sys.argv[1], 'r') as reader:
#	G = nx.DiGraph()
#	content = reader.readlines()

#	for row in content:
#		gene1, gene2, weight = row.split()
#		G.add_edge(gene1, gene2, weight=abs(float(weight)))



#take adjacency matrix instead of an edge list
with open(sys.argv[1], 'r') as f: 
	reader = csv.reader(f)
	data_as_list = list(reader)
	genes = data_as_list[0]
	del genes[0]
	del data_as_list[0]
	#print data_as_list

	for row in data_as_list:
		del row[0]
		for index in range(0, len(row)):
			row[index] = float(row[index])
	#print data_as_list


mat = np.array(data_as_list)
#print mat

G = nx.from_numpy_matrix(mat, create_using=nx.DiGraph())

#read through connected components
components = nx.weakly_connected_components(G)


#find clusters
for component in components:
	#print "curr componet", component

	#create matrix
	adjacency_matrix = nx.to_numpy_matrix(G, nodelist=component)
	#print adjacency_matrix
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
	#print "Converged to: \n", old
	#find cluster!
	#every pair (i, j) which has a non-zero value are in the same cluster
	for y in range(0, len(old)):	#col
		for z in range(0, len(old)):	#row
			if old.item(z,y) > 0.00001:
				#print "Hit", y, z
				add_flag = 0
				for c in clusterList: 
					if component[y] in c or component[z] in c:
						c.add(component[z])
						c.add(component[y])
						add_flag = 1
						break
				if add_flag == 0:
					#print "new cluster:", component[y], component[z]
					new_cluster = [component[y], component[z]]
					clusterList.append(set(new_cluster))


#print clusterList
#output textfile

writer = open('output.txt', 'w')
for cl in range(0, len(list(clusterList))):
	current_cluster = list(clusterList[cl])
	#print current_cluster
	for el in current_cluster:
		#print genes[el], cl+1
		writer.write("%s \t%s\n" % (genes[el], cl+1))



#for cl in range(1, len(clusterList)+1):
#	current_cluster = clusterList(cl-1)
#	print current_cluster
#	for el in list(current_cluster):
#		print genes[el], cl





	#calculate change of matrices
	#diff = np.subtract(new, old)
	#print diff
