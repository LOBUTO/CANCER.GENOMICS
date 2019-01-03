# cgp_encode.py
# Function for gene set deep autoencoder

import tensorflow as tf
import pandas as pd
import numpy as np
import sys
import random
import pickle
import math
from subprocess import Popen, PIPE
from sklearn import preprocessing
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
	n_th   = n.split("_")[0] # In case we have random sampling label
	counts = counts[counts.Gene_set > int(n_th)]
	g_s    = list(counts.index)

	gsea   = gsea[gsea.Gene_set.isin(g_s)]
	print("Filtered number of gene sets: %s"%(len(set(gsea.Gene_set))))

	if "random" in n:
		print("Applying randomnization to gene sets")
		all_genes   = list(set(gsea.genes))
		all_sets    = list(set(gsea.Gene_set))
		
		set_n       = {i:int(counts.loc[[i]].Gene_set) for i in all_sets}
		counts_dict = {i:random.sample(all_genes, set_n[i]) for i in set_n.keys()} #RANDOM SAMPLING
		
		counts_dict = [pd.DataFrame({"Gene_set":i, "genes":counts_dict[i]}) for i in counts_dict.keys()] #ASSEMBLE
		gsea        = pd.concat(counts_dict)

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

def train_valid_test_split(line_count):
	# Splits 80/10/10
	random.seed(1234) #NOTE Remove for total random (rnd)

	temp_index  = random.sample(xrange(line_count), line_count) #Randomize at first
	temp_index  = np.array_split(temp_index, 5)

	validation_index        = np.array_split(temp_index[4], 2)
	valid_index, test_index = list(validation_index[0]), list(validation_index[1]) #Validate on 10% and test on 10%
	train_index = list(np.concatenate(temp_index[:4])) #Train on 80%

	return train_index, valid_index, test_index

def layers(n, l):
	# n = original number of features
	# l = number of layers
    if l==1: 
        layers = [n, 1]
    else:
        layers = [n]
        for i in xrange(l-1):
            n = n/2
            layers.append(n)
            
        layers = layers + [1]
    return layers

def deep_layers(n, l):
	# n = starting features
	# l = number of layers
    
    layers = [n]
    for i in xrange(l-1):
        n = n/2
        layers.append(n)
        
    return layers

def scale_0_1_multiple(x, y, z=np.array([0])):
	# Will scale between 0-1 by column
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per gene across dataset: multiple")
	min_max_scaler = preprocessing.MinMaxScaler()
	x = min_max_scaler.fit_transform(x)
	y = min_max_scaler.transform(y)

	if (z.shape[0] > 1):
		z = min_max_scaler.transform(z)
		return x, y, z
	else:
		return x, y

def scale_0_1_single(x):
	print("Scaling per gene across dataset: single")

	min_max_scaler = preprocessing.MinMaxScaler()
	x = min_max_scaler.fit_transform(x)
	return x

def arch_layers(start, arch):
	arch   = arch.split("_")

	layers = [start]
	for i in arch:
		layers.append(start/int(i))

	return layers

def dense_batch_lrelu(x, phase, scope):
	# x: matmul(a,b), no bias needed
	with tf.variable_scope(scope):
		h = tf.contrib.layers.batch_norm(x,
										 center=True, scale=True, 
                                         is_training=phase,
                                         scope='bn')

        return lrelu(h, alpha=0.2)

# Load arguments
in_folder     = "/tigress/zamalloa/CGP_FILES/"
out_folder    = "/tigress/zamalloa/GSEA_FILES/RESULTS/"
gsea          = sys.argv[1] #c2setcover, cancer, c4_cancer
g_filter      = sys.argv[2] #50, 100, 200, 250 -or- n_random_1, n_random_2
arch          = sys.argv[3] # The architecture 2_4_8 4_8 2_4 4
keepprob      = float(sys.argv[4])

in_file       = in_folder + "070818_cgp_exp.txt"

# Load files
data_feat     = pd.read_csv(in_file, sep="\t") # samples x genes
train_index, valid_index, test_index       = train_valid_test_split(data_feat.shape[0])

train_feat    = data_feat.iloc[train_index,:]
valid_feat    = data_feat.iloc[valid_index,:]
test_feat     = data_feat.iloc[test_index,:]

gsea_table    = pd.read_csv("GSEA_FILES/{}_sets".format(gsea), sep="\t")
gsea_table    = gsea_table[gsea_table.genes.isin(list(data_feat.columns.values))]

# Prep output files
log_file      = "{}{}_{}_{}_deepautoencoder_{}_{}_log".format(out_folder, "CGP", gsea,  g_filter, arch, keepprob)
with open(log_file, "w") as f:
	f.write("Epoch\tbs\tkeep_prob\ttrain_error\tvalid_error\ttest_error")

# Parameters
slr           = 0.001 #0.0001 #0.00001
batch_size    = 20

# Execute
gsea_table    = filter_gsea(gsea_table, g_filter)
g_sets        = list(set(gsea_table.Gene_set))
all_genes     = list(set(gsea_table.genes))
train_data    = train_feat[all_genes]
valid_data    = valid_feat[all_genes]
test_data     = test_feat[all_genes]

train_feat    = np.asarray(train_data) #[samples, genes]
valid_feat	  = np.asarray(valid_data)
test_feat     = np.asarray(test_data)

# Standardize initial inputs to 0-1 range
train_feat, valid_feat, test_feat = scale_0_1_multiple(train_feat, valid_feat, test_feat)

print("SHAPES!!!")
print(train_feat.shape)
print(valid_feat.shape)
print(test_feat.shape)
print("Number of gene sets used for parallel encoding: {}".format(len(g_sets)))

# Build graphs
init_features = len(all_genes)
layers        = arch_layers(len(g_sets), arch)
rev_layers    = [i for i in reversed(layers)]
# sess          = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
main_dict     = {}

# Placeholders
with tf.device('/gpu:0'):
	x           = tf.placeholder(tf.float32, shape=[None, init_features], name="x")
	phase_train = tf.placeholder(tf.bool, name='phase_train')
	keep_prob   = tf.placeholder(tf.float32, name="keep_prob")

# ENCODING
g_loc     = {g:[train_data.columns.get_loc(i) for i in list(gsea_table[gsea_table.Gene_set==g].genes)] for g in g_sets}
g_sizes   = {g: len(set(gsea_table[gsea_table.Gene_set==g].genes)) for g in g_sets}

with tf.device('/gpu:0'):

	x_noise = tf.nn.dropout(x, keep_prob) #Apply noise to expression values
	for g in g_sets:
		n_feat                     = g_sizes[g]

		main_dict["W0_e_%s"%g] 	   = weight_variable([ n_feat ,1]) #[n_genes, half_n_genes]
		# main_dict["b0_e_%s"%g]     = bias_variable([1])

		main_dict["loc_m_%s"%g]    = tf.transpose(tf.concat([[x_noise[:,i]] for i in g_loc[g]], 0)) #[n_genes, batch]
		main_dict["h0_e_%s"%g]     = tf.matmul(main_dict["loc_m_%s"%g], main_dict["W0_e_%s"%g]) #no bias needed

	main_dict["h0_e"]  = dense_batch_lrelu(tf.concat([main_dict["h0_e_%s"%g] for g in g_sets], 1),
										   phase_train, scope="h0_e") #bn includes bias
	main_dict["h0d_e"] = tf.nn.dropout(main_dict["h0_e"], keep_prob) # Dropout

	# Additional layers after gene set layer
	l = -1 # In case there are no additional layers
	for l in xrange(len(layers)-1):
		main_dict["W%s_e"%(l+1)]   = weight_variable([ layers[l] , layers[l+1]])
		# main_dict["b%s_e"%(l+1)]   = bias_variable([layers[l+1]])

		main_dict["h%s_e"%(l+1)]   = dense_batch_lrelu(tf.matmul(main_dict["h%sd_e"%l], main_dict["W%s_e"%(l+1)]), 
														phase_train, scope="h%s_e"%(l+1)) 
		main_dict["h%sd_e"%(l+1)]  = tf.nn.dropout(main_dict["h%s_e"%(l+1)], keep_prob)

# DECODING
with tf.device('/gpu:1'):

	main_dict["h0d_d"] = main_dict["h%sd_e"%(l+1)]
		
	# Decode hidden
	for k in xrange(len(rev_layers)-1):
		main_dict["W%s_d"%(k+1)]  = weight_variable([rev_layers[k], rev_layers[k+1]])
		# main_dict["b%s_d"%(k+1)]  = bias_variable([rev_layers[k+1]])

		main_dict["h%s_d"%(k+1)]  = dense_batch_lrelu(tf.matmul(main_dict["h%sd_d"%k], main_dict["W%s_d"%(k+1)]), 
												      phase_train, scope="h%s_d"%(k+1))
		main_dict["h%sd_d"%(k+1)] = tf.nn.dropout(main_dict["h%s_d"%(k+1)], keep_prob)

	# Decode to main
	main_dict["W_d"] = weight_variable([rev_layers[-1], init_features])
	main_dict["b_d"] = bias_variable([init_features])
	main_dict["h_d"] = tf.nn.sigmoid(tf.matmul(main_dict["h%sd_d"%(k+1)], main_dict["W_d"]) + main_dict["b_d"])

# Calculate error with original expression values (denoised)
with tf.name_scope('cost'):
	cost        = tf.sqrt(tf.reduce_mean(tf.square(x - main_dict["h_d"] ))) / (tf.reduce_max(x) - tf.reduce_min(x)) # NRMSE
	coded_feat  = main_dict["h%s_e"%(l+1)]

#################################################### RUNNING ########################################################
update_ops  = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
with tf.control_dependencies(update_ops):
	# Ensures that we execute the update_ops before performing the train_step
	train_step            = tf.train.AdamOptimizer(slr).minimize(cost)

# sess = tf.Session()
sess          = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
sess.run(tf.global_variables_initializer())
# sess.run(tf.group(tf.initialize_all_variables(), tf.initialize_variables(tf.local_variables()))) 
# saver                 = tf.train.Saver(max_to_keep = None)

# Validation checks
best_valid     = np.Inf
patience       = 0 # Epoch-based
max_patience   = 50 # Epoch-based
max_training   = 8000 # Maximum number of epochs to train for
count_batch    = 0 # To get tabs of batch number
count          = 0 # Number of epochs
check          = 20 # Batch size during cost checking
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

		train_error    = [cost.eval(feed_dict =
											{x:train_feat[i*check : (i+1)*check,],
											 keep_prob:1.0, phase_train: False})
						  for i in xrange(train_checks)]
		
		train_error    = np.mean(train_error)
		# train_error    = cost.eval(feed_dict = {x:train_feat, keep_prob:1.0, phase_train:False})
		valid_loss     = cost.eval(feed_dict = {x:valid_feat, keep_prob:1.0, phase_train:False})
		test_loss      = cost.eval(feed_dict = {x:test_feat, keep_prob:1.0, phase_train:False})

		with open(log_file, "a") as f:
			f.write("\n%s\t%s\t%s\t%s\t%s\t%s"%(count, batch_size, keepprob, train_error, valid_loss, test_loss))

		# Update validation
		if valid_loss < best_valid:
			best_valid=valid_loss
			patience=0

			# UPDATE: Evaluate features for best model
			train_encoded         = coded_feat.eval(feed_dict = {x:train_feat, keep_prob:1.0, phase_train:False})
			valid_encoded         = coded_feat.eval(feed_dict = {x:valid_feat, keep_prob:1.0, phase_train:False})
			test_encoded          = coded_feat.eval(feed_dict = {x:test_feat, keep_prob:1.0, phase_train:False})
		else:
			patience+=1
	else:
		count_batch+=1
print("Done modeling CNN")

# Write out best model
train_encoded         = pd.DataFrame(train_encoded)
train_encoded.columns = [str(i) for i in xrange(layers[-1])]
train_encoded.index   = list(train_data.index)
valid_encoded         = pd.DataFrame(valid_encoded)
valid_encoded.columns = [str(i) for i in xrange(layers[-1])]
valid_encoded.index   = list(valid_data.index)
test_encoded          = pd.DataFrame(test_encoded)
test_encoded.columns  = [str(i) for i in xrange(layers[-1])]
test_encoded.index    = list(test_data.index)

train_file    = "{}{}_{}_{}_deepautoencoder_{}_{}_train.txt".format(out_folder, "CGP", gsea, g_filter, arch, keepprob)
valid_file    = "{}{}_{}_{}_deepautoencoder_{}_{}_valid.txt".format(out_folder, "CGP", gsea, g_filter, arch, keepprob)
test_file     = "{}{}_{}_{}_deepautoencoder_{}_{}_test.txt".format(out_folder, "CGP", gsea, g_filter, arch, keepprob)
train_encoded.to_csv(train_file, sep="\t", header=True, index=True)
valid_encoded.to_csv(valid_file, sep="\t", header=True, index=True)
test_encoded.to_csv(test_file, sep="\t", header=True, index=True)

sess.close()
print("DONE")