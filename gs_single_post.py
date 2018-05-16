# gs_single_post.py
# NOTE: Assumes that lowest filter setting has been ran (20) and draws all single independent autoencoders from there

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
target        = sys.argv[1]
gsea          = sys.argv[2] #cancer, c4_cancer
norm          = sys.argv[3] #sample, gene
gsea_m        = sys.argv[4] #mean
g_filter      = int(sys.argv[5]) #50, 100, 200, 250
error         = sys.argv[6] #nrmse

gee_exp       = "FALSE" # cgp GEO used by Geeleher et al. 2014 processed by us
geeproc_exp   = "TRUE" # c(T,F) # cgp GEO used by Geeleher et al. 2014 pre-processed by them
gee_target    = "TRUE" # c(T,F) # target IC50 processed by Geeleher et al. 2014
geetargetprocexp = "TRUE"

# Load files
in_file    = "{}TRAIN_MATRICES/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}"
in_file    = in_file.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m)
call_file  = "{}SINGLE_GS/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}_gfilter_{}_call"
call_file  = call_file.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m, g_filter)

# Process
gsea_table = pd.read_csv(in_file + "_gsea", sep="\t")
gsea_table = filter_gsea(gsea_table, g_filter)
g_sets     = list(set(gsea_table.Gene_set))

# g_sets     = pd.read_csv(call_file, names=["gsea"])
# g_sets     = list(g_sets.gsea)

all_train = []
all_test  = []
all_gs    = []
for g in g_sets:

	out_file = "{}{}_{}_{}_{}_{}_autoencoder_1_single_{}_{}.pickle".format(out_folder, target, gsea, norm, gsea_m, 20, error, g)

	with open(out_file, "rb") as f:
		train_encoded, test_encoded, train_index, test_index, gs_call = pickle.load(f)
		all_train.append(train_encoded)
		all_test.append(test_encoded)
		all_gs.append(gs_call)

all_train = np.concatenate(all_train, 1)
all_test  = np.concatenate(all_test, 1)

train_encoded         = pd.DataFrame(all_train)
train_encoded.columns = [all_gs]
train_encoded.index   = train_index
test_encoded          = pd.DataFrame(all_test)
test_encoded.columns  = [all_gs]
test_encoded.index    = test_index

train_file    = "{}{}_{}_{}_{}_{}_autoencoder_1_single_{}_train.txt".format(out_folder, target, gsea, norm, gsea_m, g_filter, error)
test_file     = "{}{}_{}_{}_{}_{}_autoencoder_1_single_{}_test.txt".format(out_folder, target, gsea, norm, gsea_m, g_filter, error)
train_encoded.to_csv(train_file, sep="\t", header=True, index=True)
test_encoded.to_csv(test_file, sep="\t", header=True, index=True)

print("DONE")