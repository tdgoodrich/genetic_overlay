#!/usr/bin/python

import sys
import networkx as nx
import matplotlib.pyplot as plt
import OrderedSet

def check_arcs(filename="Data/yeast_data.txt"):
    """
    Checks all arcs in yeasy_data, if a->b != b->a then we print.
    """
    infile = open(filename, "r")
    data = {}
    for line in infile:
        line = line.split()
        tail = line[0]
        head = line[2]
        score = line[4]
        # If the opposite arc is different, notify us
        if (head, tail) in data and score != data[(head,tail)]:
            print "Different scores: %s -> %s is %s, but %s -> %s is %s" % \
              (tail, head, score, head, tail, data[(head,tail)])
        else:
            data[(tail, head)] = score

def build_clean_data(filename="Data/yeast_data.txt"):
    """
    Build the network according to the original paper
    (i.e. undirected where |epsilon| > 0.08, P < 0.05).
    """
    edgelist = {}
    vertices = OrderedSet.OrderedSet()
    infile = open(filename, "r")
    for line in infile:
        line = line.split()
        tail, head, score, pvalue = (line[0], line[2], float(line[4]), float(line[6]))
        vertices.add(tail)
        vertices.add(head)
        # If a valid score
        if score != "NaN" and (score > 0.08 or score < -0.08) and pvalue < 0.05:
            # If the opposite arc has already been read
            if (head, tail) in edgelist:
                # Remove both if different positivity
                if (edgelist[(head, tail)][0] > 0) != (score > 0):
                    edgelist.pop((head, tail))
                # Else keep the edge with the better p-value
                elif (edgelist[(head, tail)][1]) > pvalue:
                    edgelist.pop((head, tail))
                    edgelist[(tail, head)] = (score, pvalue)
            # Edge we're good to add the edge
            else:
                edgelist[(tail, head)] = (score, pvalue)

    print "Vertices: %d" % len(vertices)
    print "We had %d valid edges" % len(edgelist)
    return vertices, edgelist

def build_edgelist(filename="Data/clean_yeast_data.txt"):
    edgelist = {}
    vertices = OrderedSet()
    infile = open(filename, "r")
    for line in infile:
        tail, head, score, pvalue = line.split()
        vertices.add(tail)
        vertices.add(head)
    return vertices, edgelist

def build_graph(filename="Data/clean_yeast_data.txt"):
    G = nx.Graph()
    infile = open(filename, "r")
    for line in infile:
        tail, head, score, pvalue = line.split()
        G.add_node(tail)
        G.add_node(head)
        G.add_edge(tail, head, score=score, pvalue=pvalue)
    return G

def write_edgelist(edgelist, filename="Data/clean_yeast_data.txt"):
    """
    Write the edgelist to clean_yeast_data.txt in the Code directory
    """
    outfile = open(filename, "w")
    for key in edgelist:
        outfile.write("%s %s %f %f\n" % (key[0], key[1], edgelist[key][0], edgelist[key][1]))
    outfile.close()

def write_vertices(vertices, gene_order_filename="Data/gene_order.csv"):
    gene_order_outfile = open(gene_order_filename, "w")
    for v in vertices:
        gene_order_outfile.write("%s\n" % v)
    gene_order_outfile.close()

def write_adjacency_matrices(vertices, edgelist,
  score_filename="Data/score_adjacency_matrix.csv",
  pvalue_filename="Data/pvalue_adjacency_matrix.csv"):
    score_outfile = open(score_filename, "w")
    pvalue_outfile = open(pvalue_filename, "w")
    score_row = ""
    pvalue_row = ""
    header = ",".join(["%s, " % (v) for v in vertices])
    score_outfile.write(header + "\n")
    pvalue_outfile.write(header + "\n")
    for v1 in vertices:
        score_row += ("%s, " % (v1))
        for v2 in vertices:
            score_row += " %s," % edgelist.get((v1, v2), edgelist.get((v2, v1), (0,0)))[0]
            pvalue_row += " %s," % edgelist.get((v1, v2), edgelist.get((v2, v1), (0,0)))[1]
        score_row = score_row[:-1] + "\n"
        pvalue_row = pvalue_row[:-1] + "\n"
        score_outfile.write(score_row)
        pvalue_outfile.write(pvalue_row)
        score_row = ""
        pvalue_row = ""
    score_outfile.close()
    pvalue_outfile.close()

if __name__ == "__main__":
    if sys.argv[1] == "build_clean_data":
        vertices, edgelist = build_clean_data()
        write_edgelist(edgelist)
        write_vertices(vertices)
        write_adjacency_matrices(vertices, edgelist)
    elif sys.argv[1] == "check_arcs":
        check_arcs()
    elif sys.argv[1] == "build_edgelist":
        vertices, edgelist = build_edgelist()
        write_edgelist(edgelist)
        write_vertices(vertices)
        write_adjacency_matrices(vertices, edgelist)
    elif sys.argv[1] == "plot_graph":
        vertices, edgelist = build_graph()
        plot_graph(vertices, edgelist)
    elif sys.argv[1] == "gephi":
        G = build_graph()
        nx.write_gexf(G, "Data/network.gexf")
