# gs_encode.py
# Function for modified gene set autoencoder
# Minimize one geneset at a time

import tensorflow as tf
import pandas as pd
import numpy as np
import sys
import random
import pickle
import math
from subprocess import Popen, PIPE
import cPickle as pickle

# FUNCTIONS
def weight_variable(shape):
	initial = tf.truncated_normal(shape, stddev=0.1) #Normal distributionwith std 0.1
	return tf.Variable(initial)

def bias_variable(shape):
	initial = tf.constant(0.1, shape=shape) #Uniform distribution of 0.1
	return tf.Variable(initial)

def load_gsea(g_strings, in_folder):
	# Loads gene_sets
	g_list   = g_strings.split("_")
	g_format = '{}{}_sets'
	gsea     = pd.concat([pd.read_csv(g_format.format(in_folder, g), sep="\t") for g in g_list])
	return(gsea)

def filter_gsea(gsea, n):
	print("Original number of gene sets: %s"%(len(set(gsea.Gene_set))))
	counts = pd.DataFrame(gsea.Gene_set.value_counts())
	counts = counts[counts.Gene_set > n]
	g_s    = list(counts.index)

	gsea   = gsea[gsea.Gene_set.isin(g_s)]
	print("Filtered number of gene sets: %s"%(len(set(gsea.Gene_set))))
	return(gsea)

def lrelu(x,alpha=0.1):
	return tf.maximum(alpha*x,x)

# Load arguments
in_folder     = "/tigress/zamalloa/GSEA_FILES/"
out_folder    = "/tigress/zamalloa/GSEA_FILES/RESULTS/"
target        = sys.argv[1]
gsea          = sys.argv[2] #cancer, c4_cancer
norm          = sys.argv[3] #sample, gene
gsea_m        = sys.argv[4] #mean
g_filter      = int(sys.argv[5]) #50, 100, 200, 250
error         = sys.argv[6] #nrmse
gs_call       = sys.argv[7] #name of gsea
gee_exp       = "FALSE" # cgp GEO used by Geeleher et al. 2014 processed by us
geeproc_exp   = "TRUE" # c(T,F) # cgp GEO used by Geeleher et al. 2014 pre-processed by them
gee_target    = "TRUE" # c(T,F) # target IC50 processed by Geeleher et al. 2014
geetargetprocexp = "TRUE"

in_tables     = "{}SINGLE_GS/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}_gfilter_{}_gset_{}.pickle"
in_tables     = in_tables.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m, g_filter, gs_call)

# Load files
with open(in_tables, "rb") as f:
	train_feat, test_feat, train_index, test_index = pickle.load(f)

# Prep output files
log_file      = "{}{}_{}_{}_{}_{}_autoencoder_1_single_{}_{}_log".format(out_folder,target, gsea, norm, gsea_m, g_filter, error, gs_call)
with open(log_file, "w") as f:
	f.write("Epoch\tCompound\tGene_set\tbs\tkeep_prob\ttrain_error\tvalid_error")

# Parameters
n_feat        = train_feat.shape[1]
slr           = 0.00001 #Changed from 0.0001
batch_size    = 20
keepprob      = 0.5

print("SHAPES!!!")
print(train_feat.shape)
print(test_feat.shape)

# Build graph
sess          = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
main_dict     = {}

# Placeholders
with tf.device('/gpu:0'):
	x             = tf.placeholder(tf.float32, shape=[None, n_feat], name="x")
	phase_train   = tf.placeholder(tf.bool, name='phase_train')
	keep_prob     = tf.placeholder(tf.float32, name="keep_prob")

# Encoder
with tf.device('/gpu:0'):

	main_dict["W1_e"] 	  = weight_variable([ n_feat ,1]) 
	main_dict["b1_e"]     = bias_variable([1])

	main_dict["h1_e"]     = lrelu(tf.matmul(x, main_dict["W1_e"]) + main_dict["b1_e"], alpha=0.2)

# Decoder
with tf.device('/gpu:0'):
	main_dict["W1_d"] 	  = weight_variable([1, n_feat])
	main_dict["b1_d"]     = bias_variable([n_feat])

	main_dict["y1_d"]     = tf.matmul(main_dict["h1_e"], main_dict["W1_d"]) + main_dict["b1_d"]

cost           = tf.sqrt(tf.reduce_mean(tf.square( x - main_dict["y1_d"] ))) / (tf.reduce_max(x) - tf.reduce_min(x)) # NRMSEcost
accuracy       = cost
coded_feat     = main_dict["h1_e"]

# Initialize all variables
train_step     = tf.train.AdamOptimizer(slr).minimize(cost) 
sess.run(tf.group(tf.initialize_all_variables(), tf.initialize_variables(tf.local_variables()))) 
saver          = tf.train.Saver(max_to_keep = None)

# Validation checks
check          = 10 # Batch size during accuracy checking
train_checks   = int(math.ceil(train_feat.shape[0] / float(check) )) #n_samples/check
n_batches      = train_feat.shape[0] / batch_size # n_samples/batch_size

best_valid     = np.Inf
patience       = 0 # Epoch-based
max_patience   = 30 # Epoch-based
max_training   = 3000 # Maximum number of epochs to train for
count_batch    = 0 # To get tabs of batch number
count          = 0 # Number of epochs

# Train
while (patience < max_patience) and (count < max_training):

	feat_batch   = train_feat[count_batch*batch_size : (count_batch+1)*batch_size, ]
	train_step.run(feed_dict={x: feat_batch, keep_prob: keepprob, phase_train: True})

	if count_batch == (n_batches-1): # Batch limit count, that is a whole epoch has passed
		count_batch = 0
		count+=1 # Epoch number increase

		train_error    = [accuracy.eval(feed_dict =
											{x:train_feat[i*check : (i+1)*check,],
											 keep_prob:1.0, phase_train: True})
						  for i in xrange(train_checks)]
		train_error    = np.mean(train_error)
		valid_loss     = accuracy.eval(feed_dict = {x:test_feat,
													keep_prob:1.0, phase_train:True})

		with open(log_file, "a") as f:
			f.write("\n%s\t%s\t%s\t%s\t%s\t%s"%(count, target, batch_size, keepprob, train_error, valid_loss))

		# Update validation
		if (valid_loss + 0.0001)<best_valid:
			best_valid=valid_loss
			patience=0

			# UPDATE: Evaluate features for best model
			train_encoded    = coded_feat.eval(feed_dict = {x:train_feat, keep_prob:1.0, phase_train:True})
			test_encoded     = coded_feat.eval(feed_dict = {x:test_feat, keep_prob:1.0, phase_train:True})
		else:
			patience+=1
	else:
		count_batch+=1

print("Done modeling CNN")

# Write out best model
print("Storing encoded features")
out_file = "{}{}_{}_{}_{}_{}_autoencoder_1_single_{}_{}.pickle".format(out_folder, target, gsea, norm, gsea_m, g_filter, error, gs_call)
with open(out_file, "wb") as f:
	pickle.dump([train_encoded, test_encoded, train_index, test_index, gs_call], f)

print("DONE")