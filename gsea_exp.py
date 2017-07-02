# gsea_exp.py
import sys
import pandas as pd
import numpy as np
import itertools
from scipy.stats import mannwhitneyu
from joblib import Parallel, delayed

def mann_pval(gs):

	name    = gs[0] + "$" + gs[1]
	genes_1 = gsea_dict[gs[0]]
	genes_2 = gsea_dict[gs[1]]

	genes_1_index = [i for i,n in enumerate(genes) if n in genes_1]
	genes_2_index = [i for i,n in enumerate(genes) if n in genes_2]

	gs_pvals = []
	for s in xrange(len(samples)):

		genes_1_exp = exp_matrix[genes_1_index, s]
		genes_2_exp = exp_matrix[genes_2_index, s]

		# Perform wilcoxon test
		pvals = mannwhitneyu(genes_1_exp, genes_2_exp, alternative="greater")[1]

		gs_pvals = gs_pvals + [pvals]

	# print(len(main_dict)/gene_set_length)
	return({"gs":name, "pvals":gs_pvals, "samples":samples})

# Load arguments
exp_type   = sys.argv[1]
base       = sys.argv[2] # Base can be None or cgp/ctrp
gsea_type  = sys.argv[3] #c("h", "c1", "c2.cp.biocarta", "c3", "cancer") 
both       = sys.argv[4] # T/F Both tests, not 2 sided, but one sided on both sides

# Load files
in_mac     = "/Users/jzamalloa/Documents/Rotation/PIPELINES/"
in_tiger   = "/tigress/zamalloa/"
in_lab     = "/home/zamalloa/Documents/Rotation/PIPELINES/"
in_folder  = in_tiger
gsea       = pd.read_csv(in_folder + "GSEA_FILES/" + gsea_type + "_sets", sep="\t")

if exp_type=="cgp":
	in_exp    = pd.read_csv(in_folder + "CGP_FILES/cgp_exp", sep="\t")
elif exp_type=="ctrp":
	in_exp    = pd.read_csv(in_folder + "CTRP_FILES/ctrp_exp_2.1", sep="\t")
elif exp_type=="ccle":
	in_exp    = pd.read_csv(in_folder + "CGP_FILES/ccle_exp", sep="\t")

# Pre-process
exp_matrix = in_exp.iloc[:, 1:].as_matrix()

samples    = list(in_exp)[1:]
genes      = list(in_exp["rn"])

print(gsea.shape)
gsea       = gsea[gsea.genes.isin(genes)]
gsea_dict  = {i:list(j) for i,j in gsea.groupby("Gene_set")["genes"]}
print(gsea.shape)

# Do we filter by a base?
if base=="None":
	gene_sets  = gsea.groupby(["Gene_set"]).size().reset_index(name="count")
	gene_sets  = gene_sets[gene_sets["count"] > 50]
	gene_sets  = gene_sets[gene_sets["count"] < 200]
	gene_sets  = [i for i in itertools.combinations(list(gene_sets["Gene_set"]),2)]
	
	if both=="T":
		print("both one-tailed used")
		gene_r    = [(i[1], i[0]) for i in gene_sets]
		gene_sets = gene_sets + gene_r
else:
	base       = pd.read_csv(in_folder + "GSEA_FILES/" + gsea_type + "_gsea_" + base + "_both_" + both + "_pvals", sep="\t")
	gene_sets  = list(set(base["gs"]))
	gene_sets  = [(i.split("$")[0], i.split("$")[1]) for i in gene_sets]

print(gsea_type, len(gene_sets))

# Process and write
file_out = open(in_folder + "GSEA_FILES/" + gsea_type + "_gsea_"+ exp_type + "_both_" + both + "_pvals", "w")
file_out.write("sample" + "\t" + "gs" + "\t" + "pvals")
file_out.close()

count = 0.0
gene_set_length = float(len(gene_sets))

main_dict = Parallel(n_jobs=40)(delayed(mann_pval)(i) for i in gene_sets)

print("Done calculating")
# Write to file
for gs in main_dict:
	for i in xrange(len(gs["samples"])):
		with open(in_folder + "GSEA_FILES/" + gsea_type + "_gsea_"+ exp_type + "_both_" + both + "_pvals", "a") as logfile:
			logfile.write("\n" + gs["samples"][i] + "\t" + gs["gs"] + "\t" + str(gs["pvals"][i]))


print ("Done writing")