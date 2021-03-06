import OrderedSet
import networkx as nx

################################################################################
# Misc
################################################################################

def check_arcs(filename="Data/yeast_data.txt"):
    """
    Checks all arcs in yeast_data.txt, if a->b != b->a then we print.

    Input:
        filename (string)
    """
    with open(filename, "r") as infile:
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

def convert_correlation_matrix(input_filename="Data/correlation_matrix.csv",
  output_filename="Data/Unreduced/correlation_matrix.csv"):
    """
    Converts Data/correlation_matrix.txt to Data/Unreduced/correlation_matrix.csv.
    Adds tabs, removes useless row/col.

    Input:
        filename (string)
    Output:
        Writes Data/Unreduced/correlation_matrix.csv
    """
    with open(input_filename, "r") as infile, open(output_filename, "w") as outfile:
        line = infile.readline()
        line = line.replace("\n", "")
        line = line.split("\t")
        outfile.write(("," + "".join(["%s, " % item for item in line]) + "\n"))
        infile.readline() # Throw out garbage line
        for line in infile.readlines():
            line = line.replace("\n", "").replace("NaN", "0")
            line = line.split("\t")
            del line[1]
            outfile.write(("".join(["%s, " % item for item in line]) + "\n"))

################################################################################
# Subroutines for populating data structures from files
################################################################################
def build_threshold_data(threshold, filename="Data/yeast_data.txt"):
    """
    Build the network according to the thresholds from the original paper.
    (i.e. undirected where |epsilon| > threshold, P < 0.05).

    Input:
        threshold (float)
        filename (string)
    Output:
        vertices (ordered set)
        edgelist (dictionary[(gene1, gene2)] = (score, pvalue))
    """
    edgelist = {}
    vertices = OrderedSet.OrderedSet()
    with open(filename, "r") as infile:
        print "Reading file: %s" % (filename)
        for line in infile:
            line = line.split()
            tail, head, score, pvalue = (line[0], line[2], float(line[4]), float(line[6]))
            vertices.add(tail)
            vertices.add(head)
            # If a valid score
            if score != "NaN" and (score > threshold or score < (-1 * threshold)) and pvalue < 0.05:
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

    print "Threshold %.2f gives %d nodes and %d edges" % \
      (threshold, len(vertices), len(edgelist))
    return vertices, edgelist

def build_edgelist(filename="Data/clean_yeast_data.txt"):
    """
    Read in a thresholded data file into vertices and edgelist.

    Input:
        filename (string) [not raw data]
    Output:
        vertices (ordered set)
        edgelist (dictionary[(gene1, gene2)] = (score, pvalue))
    """
    edgelist = {}
    vertices = OrderedSet.OrderedSet()
    with open(filename, "r") as infile:
        for line in infile:
            tail, head, score, pvalue = line.split()
            vertices.add(tail)
            vertices.add(head)
            edgelist[(tail, head)] = (score, pvalue)
        print "File %s has %d nodes and %d edges" % \
              (filename, len(vertices), len(edgelist))
        return vertices, edgelist

def build_graph(filename="Data/clean_yeast_data.txt"):
    """
    Read in a thresholded data file into a NetworkX graph.

    Input:
        filename (string) [not raw data]
    Output:
        G (NetworkX Graph)
    """
    G = nx.Graph()
    infile = open(filename, "r")
    for line in infile:
        tail, head, score, pvalue = line.split()
        G.add_node(tail)
        G.add_node(head)
        G.add_edge(tail, head, score=score, pvalue=pvalue)
    print "Constructed a graph with %d nodes and %d edges" % \
      (G.order(), G.size())
    return G

def build_cluster_graph(filename):
    """
    Read in a clustered data file into a NetworkX graph.

    Input:
        filename (string) [not raw data]
    Output:
        G (NetworkX Graph)
    """
    G = nx.Graph()
    infile = open(filename, "r")
    for line in infile:
        tail, head, score, pvalue, cluster = line.split()
        G.add_node(tail)
        G.add_node(head)
        G.add_edge(tail, head, score=score, pvalue=pvalue, cluster=cluster)
    print "Constructed a graph with %d nodes and %d edges" % \
      (G.order(), G.size())
    return G

################################################################################
# Subroutines for writing data structures to files
################################################################################

def write_edgelist(edgelist, filename="Data/clean_yeast_data.txt"):
    """
    Write an edgelist to a file.

    Input:
        edgelist (dictionary[(gene1, gene2)] = (score, pvalue))
        filename (string)
    Output:
        Writes the edgelist to a file with columns "gene1, gene2, score, pvalue"
    """
    with open(filename, "w") as outfile:
        print "Writing edgelist to file: %s" %(filename)
        for key in edgelist:
            outfile.write("%s %s %f %f\n" % (key[0], key[1], edgelist[key][0], edgelist[key][1]))

def write_vertices(vertices, filename="Data/gene_order.csv"):
    """
    Write the nodes to a file.

    Input:
        vertices (ordered list)
        filename (string)
    Output:
        Writes the edgelist to a file with column "gene1"
    """
    with open(filename, "w") as outfile:
        print "Writing nodes to file: %s" % (filename)
        for v in vertices:
            outfile.write("%s\n" % v)

def write_adjacency_matrices(vertices, edgelist,
  score_filename="Data/score_adjacency_matrix.csv",
  pvalue_filename="Data/pvalue_adjacency_matrix.csv"):
    """
    Writes the edgelist to adjacency matrices.

    Input:
        vertices (ordered list)
        edgelist (dictionary[(gene1, gene2)] = (score, pvalue))
        score_filename (string)
        pvalue_filename (string)
    Output:
        Write adjacency matrices to the corresponding filenames in this format:
             , gene1, gene2, gene3, ...
        gene1, score, score, score, ...
        gene2, score, score, score, ...
        ...
        One matrix uses score, the other uses pvalue.
    """

    with open(score_filename, "w") as score_outfile, open(pvalue_filename, "w") as pvalue_outfile:
        print "Writing score adjacency matrix: %s" % (score_filename)
        print "Writing pvalue adjacency matrix: %s" % (pvalue_filename)
        score_row = ""
        pvalue_row = ""
        header = "," + "".join(["%s, " % (v) for v in vertices])
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

# TODO: Hack, please replace
def write_with_clusters(vertices, edgelist, cluster_filename, output_filename):
    """
    Writes an edgelist file with clusters.
    If both endpoints are in the cluster, we know the edge's cluster.
    Otherwise we remove the edge as our reduction.

    Input:
        vertices (ordered list)
        edgelist (dictionary[(gene1, gene2)] = (score, pvalue))
        cluster_filename (string)
    Output:
        Writes the edgelist to a file with columns "gene1, gene2, score, pvalue, cluster"
    """

    cluster_lookup = {}
    cluster_data = open(cluster_filename, "r")
    for line in cluster_data.readlines():
        gene, cluster = line.split()
    #    cluster_lookup[gene] = cluster_lookup.get(gene, []) + [cluster]
        cluster_lookup[gene] = cluster
    #for key in cluster_lookup:
    #    print "Cluster ", key, " size: ", len(cluster_lookup[key])

    with open(output_filename, "w") as outfile:
        for (gene1, gene2) in edgelist:
            #print "Gene1: %s  Gene2: %s  Cluster1: %s  Cluster2: %s" % (gene1, gene2, cluster_lookup[gene1], cluster_lookup[gene2])
            if cluster_lookup[gene1] == cluster_lookup[gene2]:
                score, pvalue = edgelist[(gene1, gene2)]
                outfile.write("%s %s %s %s %s\n" % (gene1, gene2, score, pvalue, cluster_lookup[gene1][0]))
