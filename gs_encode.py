# gs_encode.py
# Function for gene set autoencoder

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

def valid_indexes(line_count):

    random.seed(1234) #NOTE Remove for total random (rnd)
    temp_index  = random.sample(xrange(line_count), line_count) #Randomize at first
    valid_index = list(np.array_split(temp_index, 10)[9]) #Validate on 10%

    train_index = list(set(temp_index) - set(valid_index)) #Train on 90%
    print("Train samples: %s, Valid samples: %s"%(len(train_index), len(valid_index)) )

    return train_index, valid_index

# Load arguments
in_folder     = "/tigress/zamalloa/GSEA_FILES/"
out_folder    = "/tigress/zamalloa/GSEA_FILES/RESULTS/"
target        = sys.argv[1]
gsea          = sys.argv[2] #cancer, c4_cancer
norm          = sys.argv[3] #sample, gene
gsea_m        = sys.argv[4] #mean
g_filter      = int(sys.argv[5]) #50, 100, 200, 250
error         = sys.argv[6] #nrmse
# g_filter      = 10
# auto_filter   = int(sys.argv[5]) #500, 200, 100, 50
gee_exp       = "FALSE" # cgp GEO used by Geeleher et al. 2014 processed by us
geeproc_exp   = "TRUE" # c(T,F) # cgp GEO used by Geeleher et al. 2014 pre-processed by them
gee_target    = "TRUE" # c(T,F) # target IC50 processed by Geeleher et al. 2014
geetargetprocexp = "TRUE"

in_file       = "{}TRAIN_MATRICES/{}_gee_exp{}_geeprocexp{}_geetarget{}_geetargetprocexp{}_gsea{}_genenorm{}_featmethod{}"
in_file       = in_file.format(in_folder, target, gee_exp, geeproc_exp, gee_target, geetargetprocexp, gsea, norm, gsea_m)

# Load files
train_feat    = pd.read_csv(in_file + "_train", sep="\t", index_col=0).transpose()
test_feat     = pd.read_csv(in_file + "_test", sep="\t", index_col=0).transpose()
gsea_table    = pd.read_csv(in_file + "_gsea", sep="\t")

# Prep output files
log_file      = "{}{}_{}_{}_{}_{}_autoencoder_1_early_{}_log".format(out_folder,target, gsea, norm, gsea_m, g_filter, error) #, auto_filter)
with open(log_file, "w") as f:
	f.write("Epoch\tCompound\tbs\tkeep_prob\ttrain_error\tvalid_error")

# Parameters
slr           = 0.00001
batch_size    = 20
keepprob      = 0.5

# Execute
gsea_table    = filter_gsea(gsea_table, g_filter)
g_sets        = list(set(gsea_table.Gene_set))
all_genes     = list(set(gsea_table.genes))
train_data    = train_feat[all_genes]
test_data     = test_feat[all_genes]

train_feat    = np.asarray(train_data) #[samples, genes]
test_feat     = np.asarray(test_data)

# Obtain validation set for early stopping
all_lines     = train_feat.shape[0]
train_index, valid_index  = valid_indexes(all_lines)
valid_feat    = train_feat[valid_index]
train_feat    = train_feat[train_index]

print(train_data.iloc[:3,:3])
print("SHAPES!!!")
print(train_feat.shape)
print(valid_feat.shape)
print(test_feat.shape)

# Build graph
sess      = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
main_dict = {}

# Placeholders
with tf.device('/gpu:0'):
	x           = tf.placeholder(tf.float32, shape=[None, len(all_genes)], name="x")
	phase_train = tf.placeholder(tf.bool, name='phase_train')
	keep_prob   = tf.placeholder(tf.float32, name="keep_prob")

# Encoder
all_genes = sum([len(set(gsea_table[gsea_table.Gene_set==g].genes)) for g in g_sets])
g_loc     = {g:[train_data.columns.get_loc(i) for i in list(gsea_table[gsea_table.Gene_set==g].genes)] for g in g_sets}

with tf.device('/gpu:1'):
	for g in g_sets:
		n_feat                     = len(set(gsea_table[gsea_table.Gene_set==g].genes))
		# half_feat                  = n_feat/2

		main_dict["W1_e_%s"%g] 	   = weight_variable([ n_feat ,1]) #[n_genes, half_n_genes]
		main_dict["b1_e_%s"%g]     = bias_variable([1])
		
		# g_loc               	   = [train_data.columns.get_loc(i) for i in list(gsea_table[gsea_table.Gene_set==g].genes)]
		main_dict["loc_m_%s"%g]    = tf.transpose(tf.concat([[x[:,i]] for i in g_loc[g]], 0)) #[n_genes, batch]

		main_dict["h1_e_%s"%g]     = lrelu(tf.matmul(main_dict["loc_m_%s"%g], main_dict["W1_e_%s"%g]) + main_dict["b1_e_%s"%g], alpha=0.2)
		# main_dict["h1d_e_%s"%g]    = tf.nn.dropout(main_dict["h1_e_%s"%g], keep_prob)

		# main_dict["W2_e_%s"%g] 	   = weight_variable([ half_feat , 1]) #[half_n_genes, 1]
		# main_dict["b2_e_%s"%g]     = bias_variable([1])
		# main_dict["h2_e_%s"%g]     = lrelu(tf.matmul(main_dict["h1d_e_%s"%g], main_dict["W2_e_%s"%g]) + main_dict["b2_e_%s"%g], alpha=0.2)
	# main_dict["W_concat"]   = tf.concat([main_dict["h1_e_%s"%g] for g in g_sets], 1)
	
	# main_dict["W2_e"]       = weight_variable([ len(g_sets) , auto_filter])
	# main_dict["b2_e"]       = bias_variable([auto_filter])
	# main_dict["h2_e"]       = lrelu(tf.matmul(main_dict["W_concat"], main_dict["W2_e"]) + main_dict["b2_e"], alpha=0.2) #LATENT!!!!
	# main_dict["h2d_e"]      = tf.nn.dropout(main_dict["h2_e"], keep_prob)

# Decoder
with tf.device('/gpu:1'):
	# main_dict["W1_d"]       = weight_variable([auto_filter, len(g_sets)])
	# main_dict["b1_d"]       = bias_variable([len(g_sets)])
	# main_dict["y1_d"]       = lrelu(tf.matmul(main_dict["h2d_e"], main_dict["W1_d"]) + main_dict["b1_d"], alpha=0.2) #All g_sets size
	# main_dict["y1d_d"]      = tf.nn.dropout(main_dict["y1_d"], keep_prob)

	# main_dict["W2_d"]       = weight_variable([len(g_sets), all_genes])
	# main_dict["b2_d"]       = bias_variable([all_genes])
	# main_dict["y2_d"]       = lrelu(tf.matmul(main_dict["y1d_d"], main_dict["W2_d"]) + main_dict["b2_d"], alpha=0.2) #All genes size

	for g in g_sets:
		n_feat                     = len(set(gsea_table[gsea_table.Gene_set==g].genes))
	# 	half_feat                  = n_feat/2

		main_dict["W1_d_%s"%g] 	   = weight_variable([1, n_feat]) # [1, n_genes]
		main_dict["b1_d_%s"%g]     = bias_variable([n_feat]) # [n_genes]

		main_dict["y1_d_%s"%g]     = tf.matmul(main_dict["h1_e_%s"%g], main_dict["W1_d_%s"%g]) + main_dict["b1_d_%s"%g] #[n_genes]
	# 	main_dict["y1d_d_%s"%g]    = tf.nn.dropout(main_dict["y1_d_%s"%g], keep_prob)

	# 	main_dict["W2_d_%s"%g] 	   = weight_variable([half_feat, n_feat]) # [1, n_genes]
	# 	main_dict["b2_d_%s"%g]     = bias_variable([n_feat]) # [n_genes]
	# 	main_dict["y2_d_%s"%g]     = tf.matmul(main_dict["y1d_d_%s"%g], main_dict["W2_d_%s"%g]) + main_dict["b2_d_%s"%g]

# Calculate error
for g in g_sets:
	main_dict["loc_m_%s"%g]  = tf.transpose(tf.concat([[x[:,i]] for i in g_loc[g]], 0))
	# main_dict["y2_d_%s"%g]   = tf.transpose(tf.concat([[main_dict["y2_d"][:,i]] for i in g_loc], 0)) #NEW
	main_dict["cost_%s"%g]   = tf.sqrt(tf.reduce_mean(tf.square( main_dict["loc_m_%s"%g] - main_dict["y1_d_%s"%g] ))) / (tf.reduce_max(main_dict["loc_m_%s"%g]) - tf.reduce_min(main_dict["loc_m_%s"%g])) # NRMSE
	

cost        = tf.reduce_mean([main_dict["cost_%s"%g] for c in g_sets])
accuracy    = cost
coded_feat  = tf.concat([main_dict["h1_e_%s"%g] for g in g_sets], 1)
# coded_feat  = main_dict["h2_e"]

# Run
train_step            = tf.train.AdamOptimizer(slr).minimize(cost)

sess.run(tf.group(tf.initialize_all_variables(), tf.initialize_variables(tf.local_variables()))) 
saver                 = tf.train.Saver(max_to_keep = None)

# Validation checks
best_valid     = np.Inf
patience       = 0 # Epoch-based
max_patience   = 30 # Epoch-based
max_training   = 4000 # Maximum number of epochs to train for
count_batch    = 0 # To get tabs of batch number
count          = 0 # Number of epochs
check          = 10 # Batch size during accuracy checking
train_checks   = int(math.ceil(train_feat.shape[0] / float(check) )) #n_samples/check
n_batches      = train_data.shape[0] / batch_size # n_samples/batch_size

# Train
print("n_batches: %s"%n_batches)
print("train_checks: %s"%train_checks)
while (patience < max_patience) and (count < max_training):

	feat_batch   = train_feat[count_batch*batch_size : (count_batch+1)*batch_size, ]
	train_step.run(feed_dict={x: feat_batch, keep_prob: keepprob, phase_train: True})

	if count_batch == (n_batches-1): # Batch limit count, that is a whole epoch has passed
		count_batch = 0
		count+=1 # Epoch number increase
		print(count)

		train_error    = [accuracy.eval(feed_dict =
											{x:train_feat[i*check : (i+1)*check,],
											 keep_prob:1.0, phase_train: True})
						  for i in xrange(train_checks)]
		print(train_error)
		train_error    = np.mean(train_error)
		print(train_error)
		valid_loss     = accuracy.eval(feed_dict = {x:valid_feat,
													keep_prob:1.0, phase_train:True})
		print(valid_loss)
		with open(log_file, "a") as f:
			f.write("\n%s\t%s\t%s\t%s\t%s\t%s"%(count, target, batch_size, keepprob, train_error, valid_loss))

		# Update validation
		if (valid_loss + 0.0001)<best_valid:
			best_valid=valid_loss
			patience=0

			# UPDATE: Evaluate features for best model
			train_encoded         = coded_feat.eval(feed_dict = {x:np.concatenate((train_feat, valid_feat), axis=0), 
																	keep_prob:1.0, phase_train:False})
			test_encoded          = coded_feat.eval(feed_dict = {x:test_feat, keep_prob:1.0, phase_train:False})
		else:
			patience+=1
	else:
		count_batch+=1
print("Done modeling CNN")

# Write out best model
print("Storing encoded features")
train_encoded         = pd.DataFrame(train_encoded)
train_encoded.columns = [str(i) for i in xrange(len(g_sets))] #len(g_sets)
train_encoded.index   = list(train_data.index)
test_encoded          = pd.DataFrame(test_encoded)
test_encoded.columns  = [str(i) for i in xrange(len(g_sets))] #len(g_sets), auto_filter
test_encoded.index    = list(test_data.index)

# train_file    = "{}{}_{}_{}_{}_{}_autoencoder_2_train.pickle".format(out_folder,target, gsea, norm, gsea_m, g_filter)
# test_file     = "{}{}_{}_{}_{}_{}_autoencoder_2_test.pickle".format(out_folder,target, gsea, norm, gsea_m, g_filter)

# with open(train_file, "wb") as f:
# 	pickle.dump(train_encoded, f)
# with open(test_file, "wb") as f:
# 	pickle.dump(test_encoded, f)

train_file    = "{}{}_{}_{}_{}_{}_autoencoder_1_early_{}_train.txt".format(out_folder,target, gsea, norm, gsea_m, g_filter, error) #, auto_filter)
test_file     = "{}{}_{}_{}_{}_{}_autoencoder_1_early_{}_test.txt".format(out_folder,target, gsea, norm, gsea_m, g_filter, error) #, auto_filter)
train_encoded.to_csv(train_file, sep="\t", header=True, index=True)
test_encoded.to_csv(test_file, sep="\t", header=True, index=True)

print("DONE")