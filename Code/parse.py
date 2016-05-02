def check_arcs():
    """
    Checks all arcs in yeasy_data, if a->b != b->a then we print.
    """
    infile = open("Code/yeast_data.txt", "r")
    dictionary = {}
    for line in infile:
        line = line.split()
        tail = line[0]
        head = line[2]
        score = line[4]
        # If the opposite arc is different, notify us
        if (head, tail) in dictionary and score != dictionary[(head,tail)]:
            print "Different scores: %s -> %s is %s, but %s -> %s is %s" % \
              (tail, head, score, head, tail, dictionary[(head,tail)])
        else:
            dictionary[(tail, head)] = score

if __name__ == "__main__":
    check_arcs()
