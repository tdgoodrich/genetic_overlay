import go_enrich
import numpy as np

import multiprocessing as mp
from multiprocessing import Pool

import random

go_dict = go_enrich.GoDictionary()


def common_enrichment(gene1, gene2):
    set_gene1 = set(go_dict.get_terms_for_gene(gene_name=gene1))
    set_gene2 = set(go_dict.get_terms_for_gene(gene_name=gene2))
    return min(len(set_gene1.intersection(set_gene2)), 1)


def calc_score(gene_obj):
    i, j = gene_obj
    return common_enrichment(gene_names[int(i)], gene_names[int(j)])


class Evaluate():
    def pairwise_network_enrichment(self, network_matrix):
        # compute pairwise network enrichment
        edge_list, no_edge_list = self.get_edge_lists(network_matrix)

        pool = mp.Pool()
        edge_score = list(pool.imap_unordered(calc_score, edge_list))
        pool.close()
        pool.join()

        pool = mp.Pool()
        no_edge_score = list(pool.imap_unordered(calc_score,
                                                 random.sample(no_edge_list, len(edge_list))))
        pool.close()
        pool.join()

        return (len(edge_score), sum(edge_score),
                len(no_edge_score), sum(no_edge_score))

    def get_edge_lists(self, network_matrix):
        edge_list = []
        no_edge_list = []
        for i in range(network_matrix.shape[0]-1):
            for j in range(i+1, network_matrix.shape[1]):
                if i % 1000 == 1 and j == i+1:
                    print(i)

                if network_matrix[i, j] == 0:
                    no_edge_list.append((i, j))

                else:
                    edge_list.append((i, j))

        return edge_list, no_edge_list



if __name__ == "__main__":
    reduced = ['YML106W', 'YKL135C', 'YDR516C', 'YLR420W',
               'YNL111C', 'YHR007C', 'YLR014C', 'YKL216W',
               'YNL078W', 'YJR005W', 'YJL130C']

    adj_matrix = np.genfromtxt('../Data/score_adjacency_matrix.csv',
                               delimiter=',', skip_header=1)
    adj_matrix = np.delete(adj_matrix, 0, 1)

    with open('../Data/score_adjacency_matrix.csv') as f:
        gene_names = f.readline().rstrip().replace(' ', '').split(',')[1:-1]

    evaluator = Evaluate()
    print(evaluator.pairwise_network_enrichment(adj_matrix))
    # do reduced matrix?
