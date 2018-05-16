# gs_single_prep.py
# Prep files for gs_encode_single.py

import pandas as pd
import numpy as np
import sys
import random
import pickle
import math
from subprocess import Popen, PIPE
import cPickle as pickle

def filter_gsea(gsea, n):
	print("Original number of gene sets: %s"%(len(set(gsea.Gene_set))))
	counts = pd.DataFrame(gsea.Gene_set.value_counts())
	counts = counts[counts.Gene_set > n]
	g_s    = list(counts.index)

	gsea   = gsea[gsea.Gene_set.isin(g_s)]

	# Thresholding names gene_set to be able to store them as files
	gsea.Gene_set = [i[:100] for i in gsea.Gene_set]
	print("Filtered number of gene sets: %s"%(len(set(gsea.Gene_set))))
	return(gsea)

# Load arguments
in_folder     = "/tigress/zamalloa/GSEA_FILES/"
out_folder    = "/tigress/zamalloa/GSEA_FILES/RESULTS/"
target        = sys.argv[1] #docetaxel
gsea          = sys.argv[2] #cancer, c4_cancer
norm          = sys.argv[3] #sample, gene
gsea_m        = sys.argv[4] #mean
g_filter      = int(sys.argv[5]) #50, 100, 200, 250
error         = sys.argv[6] #nrmse

gee_exp       = "FALSE" # cgp GEO used by Geeleher et al. 2014 processed by us
geeproc_exp   = "TRUE" # c(T,F) # cgp GEO used by Geeleher et al. 2014 pre-processed by them
gee_target    = "TRUE" # c(T,F) # target IC50 processed by Geeleher et al. 2014
geetargetprocexp = "TRUE"

in_file       = "{}TRAIN_MATRICES/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}"
in_file       = in_file.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m)

# Load files
train_feat    = pd.read_csv(in_file + "_train", sep="\t", index_col=0).transpose()
test_feat     = pd.read_csv(in_file + "_test", sep="\t", index_col=0).transpose()
train_index   = list(train_feat.index)
test_index    = list(test_feat.index)
gsea_table    = pd.read_csv(in_file + "_gsea", sep="\t")

# Build general tables
gsea_table    = filter_gsea(gsea_table, g_filter)
g_sets        = list(set(gsea_table.Gene_set))

# Write out table per gene set
file_out = "{}SINGLE_GS/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}_gfilter_{}"
file_out = file_out.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m, g_filter)

print("Writing out single gset train and test tables:")
for g in g_sets:
	print(g)
	genes       = list(set(gsea_table[gsea_table.Gene_set==g].genes))
	train_data  = np.asarray(train_feat[genes])
	test_data   = np.asarray(test_feat[genes])

	pickle_file = file_out + "_gset_" + g + ".pickle"

	with open(pickle_file, "wb") as f:
		pickle.dump([train_data, test_data, train_index, test_index], f)

print("Writing gsea call list")
gsets = pd.DataFrame(g_sets)
gsets.to_csv(file_out + "_call", index=False, header=False)

print("Done writing single tables")