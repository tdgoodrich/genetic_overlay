import csv
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations

from collections import defaultdict
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj


class GoEnrich():
    def __init__(self):
        obodag = GODag("../Data/evaluation_reference/goslim_yeast.obo")
        background = [line.strip() for line in open('../Data/evaluation_reference/gene_list.txt')]
        geneid2gos_yeast = read_associations('../Data/evaluation_reference/geneid2gos_yeast.txt')

        self.goeaobj = GOEnrichmentStudy(
            background,
            geneid2gos_yeast,
            obodag,
            propogate_counts=False,
            alpha=0.05,
            methods=['fdr_bh'])

    def measure_enrichment(self,
                           gene_set=['YML106W', 'YKL135C', 'YDR516C',
                                     'YLR420W', 'YNL111C', 'YHR007C',
                                     'YLR014C', 'YKL216W', 'YNL078W',
                                     'YJR005W', 'YJL130C'],
                           run_name='base',
                           cluster_id=1):

        gene_ids = ['YML106W', 'YKL135C', 'YDR516C', 'YLR420W', 'YNL111C',
                    'YHR007C', 'YLR014C', 'YKL216W', 'YNL078W', 'YJR005W',
                    'YJL130C']

        goea_results_all = self.goeaobj.run_study(gene_ids)

        # we can get significant only
        # goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

        self.goeaobj.wr_txt("../Results/" + run_name + "_" + str(cluster_id) +
                            ".txt", goea_results_all)


class GoDictionary():
    def __init__(self):
        self.go_dict = self.parse_terms()
        # functional -
        # self.go_dict = self.parse_functional_terms()

    def get_terms_for_gene(self, gene_name='YBL021C'):
        if gene_name in self.go_dict.keys():
            return self.go_dict[gene_name]
        else:
            return []

    def gene_association_list(self, filename='../Data/evaluation_reference/go_slim_mapping.tab'):
        gene_list = defaultdict(set)

        output_file = open('../Data/evaluation_reference/geneid2gos_yeast.txt', 'w')

        with open(filename, 'r') as csvfile:
            genereader = csv.reader(csvfile, delimiter='\t')
            for row in genereader:
                if row[3] == 'P' and row[5]:
                    output_file.write(row[0] + '\t' + row[5] + '\n')

    # returning only GO id for now
    def parse_terms(self, filename='../Data/evaluation_reference/go_slim_mapping.tab'):
        gene_dict = {}
        with open(filename, 'r') as csvfile:
            genereader = csv.reader(csvfile, delimiter='\t')
            for row in genereader:
                if row[3] == 'P':
                    if row[0] in gene_dict.keys():
                        gene_dict[row[0]].append(row[5])
                    else:
                        gene_dict[row[0]] = [row[5]]
        return gene_dict

    def parse_functional_terms(self,
                               go_mapping='../Data/evaluation_reference/go_slim_mapping.tab',
                               go_functional_filename='../Data/evaluation_reference/GO_functional_slim.txt'):

        go_functional_list = []
        with open(go_functional_filename, 'r') as csvfile:
            goreader = csv.reader(csvfile, delimiter='\t')
            for row in goreader:
                go_functional_list.append(row[0])

        gene_dict = self.parse_terms(go_mapping)

        functional_dict = {}
        for key in gene_dict:
            functional_dict[key] = []
            for annotation in gene_dict[key]:
                if annotation in go_functional_list:
                    functional_dict[key].append(annotation)

        return functional_dict

if __name__ == "__main__":
    # go_dict = GoDictionary()
    # go_dict.gene_association_list()

    # print(go_dict.get_terms_for_gene('YML106W'))
    # print(go_dict.get_terms_for_gene('YKL135C'))
    # print(go_dict.get_terms_for_gene('YDR516C'))
    # print(go_dict.get_terms_for_gene('YLR420W'))
    # print(go_dict.get_terms_for_gene('YNL111C'))
    # print(go_dict.get_terms_for_gene('YHR007C'))
    # print(go_dict.get_terms_for_gene('YLR014C'))
    # print(go_dict.get_terms_for_gene('YKL216W'))
    # print(go_dict.get_terms_for_gene('YNL078W'))
    # print(go_dict.get_terms_for_gene('YJL130C'))

    go_enrich = GoEnrich()
    go_enrich.measure_enrichment()
