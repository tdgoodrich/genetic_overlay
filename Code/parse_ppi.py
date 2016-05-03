import sys

YEAST_TAX_ID = '559292'


def build_graph(filename):
    """
    Build the network according to the original paper
    (i.e. undirected where |epsilon| > 0.08, P < 0.05).
    """
    edgelist = {}
    vertices = set()
    infile = open(filename, "r")
    for line in infile:
        line = line.split('\t')
        head, tail, head_tax, tail_tax = (line[5], line[6], line[15], line[16])
        if head_tax != YEAST_TAX_ID or tail_tax != YEAST_TAX_ID:
            continue
        vertices.add(tail)
        vertices.add(head)
        # symmetric, dumping both into the sparse dict at the moment
        edgelist[(tail, head)] = 1
        edgelist[(head, tail)] = 1

    infile.close()
    print "Vertices: %d" % len(vertices)
    print "We had %d valid edges" % len(edgelist)
    return vertices, edgelist


def write_edgelist(edgelist, filename):
    """
    Write a sparse format list of edges to sparse-yeast-ppi.txt.

    Uses edgelist - if symmetric edgelist should specifically include both
    directions.
    """
    outfile = open(filename, "w")
    for key in edgelist:
        outfile.write("%s %s %i\n" % (key[0], key[1], edgelist[key]))
    outfile.close()


def write_adjacency_matrices(edgelist, filename, vertex_file):
    """
    Write a dense format matrix of edges to dense-yeast-ppi.txt.

    Uses edgelist - if symmetric edgelist should specifically include both
    directions. If passed, it will align output to a file of vertices. If an
    identifier doesn't exist in that file but does in edgelist, it will be
    skipped. If it does exist in that file and doesn't in edgelist, all edges
    will be 0.
    """
    dense_outfile = open(filename, "w")
    dense_row = ""
    for v1 in vertices:
        for v2 in vertices:
            try:
                val = edgelist[(v1, v2)]
            except KeyError:
                val = 0
            dense_row += "%s," % val
        dense_row = dense_row[:-1] + "\n"
        dense_outfile.write(dense_row)
        dense_row = ""
    dense_outfile.close()

if __name__ == "__main__":
    if sys.argv[1] == "build_graph":
        vertices, edgelist = build_graph("Data/BIOGRID-MV-Physical-3.4.136.tab2.txt")
        write_edgelist(edgelist, "Data/sparse-yeast-ppi.txt")
        write_adjacency_matrices(edgelist,
                                 "Data/dense-yeast-ppi.txt",
                                 "Data/gene_order.csv")
