# gsea_filter.py
import pandas as pd

x = pd.read_csv("/Users/jzamalloa/Documents/Rotation/PIPELINES/GSEA_FILES/c2_sets", sep="\t")
y = list(x["genes"].groupby(x["Gene_set"]))

file_out = open("/Users/jzamalloa/Documents/Rotation/PIPELINES/GSEA_FILES/c2_sets_commonality", "w")
file_out.write("gsea_1" + "\t" + "gsea_2" + "\t" + "len_1" + "\t" + "len_2" + "\t" + "Score")

# Find similarity across lists
for i in y:
	print(i[0])
	for j in y:

		name_1, name_2   = i[0], j[0]
		genes_1, genes_2 = list(i[1]), list(j[1])
		len_1, len_2     = len(genes_1), len(genes_2)

		commonality      = len(set(genes_1) & set(genes_2)) / float(len(set(genes_1) | set(genes_2)))

		file_out.write("\n" + name_1 + "\t" + name_2 + "\t" + str(len_1) + "\t" + str(len_2) + "\t" + str(commonality))

file_out.close()