import go_enrich
import ppi_enrich
import numpy as np

import multiprocessing as mp
from multiprocessing import Pool

import random
import sys

go_dict = go_enrich.GoDictionary()
ppi_dict = ppi_enrich.PPIDictionary().parse_ppi()

with open('../Data/score_adjacency_matrix.csv') as f:
    gene_names = f.readline().rstrip().replace(' ', '').split(',')[1:-1]

def common_enrichment(gene1, gene2):
    set_gene1 = set(go_dict.get_terms_for_gene(gene_name=gene1))
    set_gene2 = set(go_dict.get_terms_for_gene(gene_name=gene2))
    return min(len(set_gene1.intersection(set_gene2)), 1)


# def common_ppi(gene1, gene2):
def calc_score(gene_obj):
    i, j = gene_obj
    return common_enrichment(gene_names[int(i)], gene_names[int(j)])


def calc_ppi(gene_obj):
    i, j = gene_obj
    return ppi_enrichment(i, j)


def ppi_enrichment(gene1, gene2):
    gene1_name = gene_names[gene1]
    gene2_name = gene_names[gene2]

    if gene1_name in ppi_dict.keys():
        if gene2_name in ppi_dict[gene1_name]:
            return 1
        else:
            return 0
    elif gene2_name in ppi_dict.keys():
        if gene1_name in ppi_dict[gene2_name]:
            return 1
        else:
            return 0
    else:
        return 0


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

        return len(edge_score), sum(edge_score), sum(no_edge_score)

    def ppi_network_enrichment(self, network_matrix, gene_names):
        edge_list, no_edge_list = self.get_edge_lists(network_matrix)

        pool = mp.Pool()
        edge_score = list(pool.imap_unordered(calc_ppi, edge_list))
        pool.close()
        pool.join()

        pool = mp.Pool()
        no_edge_score = list(pool.imap_unordered(calc_ppi, random.sample(no_edge_list, len(edge_list))))
        pool.close()
        pool.join()

        return len(edge_score), sum(edge_score), sum(no_edge_score)


    def get_edge_lists(self, network_matrix):
        edge_list = []
        no_edge_list = []
        for i in range(network_matrix.shape[0]-1):
            for j in range(i+1, network_matrix.shape[1]):
                # if i % 1000 == 1 and j == i+1:
                #     print(i)

                if network_matrix[i, j] == 0:
                    no_edge_list.append((i, j))

                else:
                    edge_list.append((i, j))

        return edge_list, no_edge_list



if __name__ == "__main__":
    reduced = ['YML106W', 'YKL135C', 'YDR516C', 'YLR420W',
               'YNL111C', 'YHR007C', 'YLR014C', 'YKL216W',
               'YNL078W', 'YJR005W', 'YJL130C']

    pos_neg = 0 # 0 for pos
    if sys.argv[1] == 'pos':
        pos_neg = 0
    elif sys.argv[1] == 'neg':
        pos_neg = 1

    if pos_neg == 0:
        threshold = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                     0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    else:
        threshold = [0, -0.01, -0.02, -0.03, -0.04, -0.05, -0.06, -0.07, -0.08,
                     -0.09, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8,
                     -0.9, -1]

    go_score = []
    ppi_score = []

    for t in threshold:
        adj_matrix = np.genfromtxt('../Data/score_adjacency_matrix.csv',
                                   delimiter=',', skip_header=1)
        adj_matrix = np.delete(adj_matrix, 0, 1)

        if pos_neg == 0:
            below_threshold_indices = adj_matrix < t
        # negative
        if pos_neg == 1:
            below_threshold_indices = adj_matrix > t
        adj_matrix[below_threshold_indices] = 0

        if pos_neg == 0:
            zero = adj_matrix < 0
        if pos_neg == 1:
            zero = adj_matrix > 0

        adj_matrix[zero] = 0

        evaluator = Evaluate()
        print('Threshold', t)
        new_go_score = evaluator.pairwise_network_enrichment(adj_matrix)
        new_ppi_score = evaluator.ppi_network_enrichment(adj_matrix, gene_names)

        go_score.append(new_go_score)
        ppi_score.append(new_ppi_score)

        print(new_go_score)
        print(new_ppi_score)

    print(go_score)
    print(ppi_score)
    # do reduced matrix?
