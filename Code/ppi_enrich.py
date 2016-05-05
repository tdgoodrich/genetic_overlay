import csv


class PPIDictionary():
    def parse_ppi(self, filename='../Data/sparse-yeast-ppi.txt'):
        protein_dict = {}
        with open(filename, 'r') as csvfile:
            protein_reader = csv.reader(csvfile, delimiter='\t')

            for row in protein_reader:
                row = row[0].split(' ')
                if row[0] in protein_dict.keys():
                    protein_dict[row[0]].append(row[1])
                else:
                    protein_dict[row[0]] = [row[1]]
        return protein_dict

if __name__ == "__main__":
    ppi_dict = PPIDictionary()
    protein_dict = ppi_dict.parse_ppi()
    print(len(protein_dict))
