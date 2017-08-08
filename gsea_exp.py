# gsea_exp.py
# cgp_site_gmv.R was used to build cgp_site_gmv_exp
# tcga expression sets (including tcga_site_gmv were constructed using tcga_exp.R)

import sys
import pandas as pd
import numpy as np
import itertools
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
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
		# pvals = ttest_ind(genes_1_exp, genes_2_exp, equal_var=False)[1]

		gs_pvals = gs_pvals + [pvals]

	return({"sample":samples, "gs":name, "pvals":gs_pvals})

# Load arguments
exp_type   = sys.argv[1]
base       = sys.argv[2] # Base can be None or cgp/ctrp
gsea_type  = sys.argv[3] # ("h", "c1", "c2.cp.biocarta", "c3", "cancer") - OR - joined combination (i.e. h_c1)
both       = sys.argv[4] # T/F Both tests, not 2 sided, but one sided on both sides
ext_gmv    = sys.argv[5] # T/F apply gmv normalization in this script
drug       = sys.argv[6] # Cisplatin, Bortezomib (unlike bortezomib_a in 'exp_type')

# Load files
in_mac     = "/Users/jzamalloa/Documents/Rotation/PIPELINES/"
in_tiger   = "/tigress/zamalloa/"
in_lab     = "/home/zamalloa/Documents/Rotation/PIPELINES/"
in_folder  = in_lab

gsea       = pd.DataFrame()
for g in gsea_type.split("_"):
	gsea       = pd.concat([gsea,
							pd.read_csv(in_folder + "GSEA_FILES/" + g + "_sets", sep="\t")]) 

if exp_type=="ctrp":
	in_exp    = pd.read_csv(in_folder + "CTRP_FILES/ctrp_exp_2.1", sep="\t")
else:
	in_exp    = pd.read_csv(in_folder + "CGP_FILES/" + exp_type + "_exp", sep="\t")	

	# Selecting relevant samples
	if "cgp" not in exp_type:
		target    = pd.read_csv(in_folder + "CGP_FILES/" + exp_type + "_target", sep="\t")
		target    = target[target.Compound==drug]
		in_exp    = in_exp.loc[:, ["rn"]+list(target["cell_name"])]

# Pre-process
exp_matrix = in_exp.iloc[:, 1:].as_matrix()
exp_mean   = np.mean(exp_matrix, axis=1).reshape(exp_matrix.shape[0], 1)
exp_var    = np.var(exp_matrix, axis=1, ddof=1).reshape(exp_matrix.shape[0], 1)

if ext_gmv == "T":

	print("Applying gmv normalization")
	exp_matrix = (exp_matrix - exp_mean) / exp_var

samples    = list(in_exp)[1:]
genes      = list(in_exp["rn"])

print(gsea.shape)
gsea       = gsea[gsea.genes.isin(genes)]
gsea_dict  = {i:list(j) for i,j in gsea.groupby("Gene_set")["genes"]}
print(gsea.shape)

# Do we filter by a base?
if base=="None":
	gene_sets  = gsea.groupby(["Gene_set"]).size().reset_index(name="count")
	print(gene_sets)
	gene_sets  = gene_sets[gene_sets["count"] >= 100] #Changed back to 100!!!
	print(gene_sets)
	# gene_sets  = gene_sets[gene_sets["count"] <= 150]
	# print(gene_sets)
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

main_dict = Parallel(n_jobs=40)(delayed(mann_pval)(i) for i in gene_sets)

print("Done calculating")
# Write to file
main_dict = pd.concat([pd.DataFrame(i) for i in main_dict])
file_out  = in_folder + "GSEA_FILES/" + gsea_type + "_gsea_"+ exp_type + "_both_" + both + "_ext_gmv_" + ext_gmv + "_pvals"

main_dict.to_csv(file_out, sep="\t", header=True, index=False)

print ("Done writing")