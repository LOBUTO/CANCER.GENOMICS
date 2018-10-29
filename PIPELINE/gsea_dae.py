# Set of Functions for Gene-set autoencoder
def load_gsea_updated(gsea_choice, in_folder):
	g_set = gsea_choice.split("_")

	all_g = []
	for g in g_set:
		file_in = "{}GSEA_FILES/{}_sets".format(in_folder, g)
		all_g.append(pd.read_csv(file_in, sep="\t"))

	all_g = pd.concat(all_g)

def gene_filter(gsea, train_exp, test_exp, filter):
	# Filters gsea and expression matrix ([samples x genes])

	# Find common genes
	common_genes   = set(train_exp.columns).intersection(list(test_exp.columns))
	common_genes   = list(set(gsea.genes).intersection(common_genes))

	# Apply threshold
	gsea           = gsea.loc[gsea.genes.isin(common_genes)]
	filtered_gsea  = (gsea.Gene_set.value_counts()
					  .reset_index(name="count")
					  .query("count>=@filter")
					  .rename(columns={"index":"Gene_set"}))
	filtered_gsea  = gsea.loc[gsea.Gene_set.isin(list(filtered_gsea.Gene_set))]

	# Clean up
	filtered_genes = list(set(filtered_gsea.genes))
	train_exp      = train_exp[filtered_genes]
	test_exp       = test_exp[filtered_genes]

	# Return
	print("Total feature space: {}".format(train_exp.shape[1]))
	return gsea, train_exp, test_exp

def dae(cell_train, cell_test, gsea, arch):
	# Trains autoencoder using tensorflow 

	import tensorflow as tf
	import pandas as pd
	import numpy as np
	import random, math
	from sklearn import preprocessing

	# Internal functions
	def scale_0_1_multiple(x, y, z=np.array([0])):
		# Will scale between 0-1 by column
		# Scaling parameters will be obtained for train(x), 
		#  and then applied to train and test data(y)

		print("Scaling per gene across dataset: multiple")
		min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
		x = min_max_scaler.fit_transform(x)
		y = min_max_scaler.transform(y)

		if (z.shape[0] > 1):
			z = min_max_scaler.transform(z)
			return x, y, z
		else:
			return x, y

	def weight_variable(shape):
		initial = tf.truncated_normal(shape, stddev=0.1) #Normal distributionwith std 0.1
		return tf.Variable(initial)

	def bias_variable(shape):
		initial = tf.constant(0.1, shape=shape) #Uniform distribution of 0.1
		return tf.Variable(initial)

	def lrelu(x,alpha=0.1):
		return tf.maximum(alpha*x,x)

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

	# Hyperparameters
	slr           = 0.0001 #0.00001
	batch_size    = 20
	keepprob      = 0.5

	# Data
	g_sets        = list(set(gsea.Gene_set))
	all_genes     = list(set(gsea.genes))
	train_feat    = np.asarray(cell_train) #[samples, genes]
	test_feat     = np.asarray(cell_test)
	log_table     = []

	# Standardize initial inputs to 0-1 range
	train_feat, test_feat = scale_0_1_multiple(train_feat, test_feat)

	print("Number of gene sets used for parallel encoding: {}".format(len(g_sets)))

	# Build Tensorflow Graphs
	init_features = len(all_genes)
	layers        = arch_layers(len(g_sets), arch)
	rev_layers    = [i for i in reversed(layers)]
	main_dict     = {}

	# Placeholders
	with tf.device('/gpu:0'):
		x           = tf.placeholder(tf.float32, shape=[None, init_features], name="x")
		phase_train = tf.placeholder(tf.bool, name='phase_train')
		keep_prob   = tf.placeholder(tf.float32, name="keep_prob")

	# ENCODING
	g_loc     = {g:[cell_train.columns.get_loc(i) for i in list(gsea[gsea.Gene_set==g].genes)] for g in g_sets}
	g_sizes   = {g: len(set(gsea[gsea.Gene_set==g].genes)) for g in g_sets}

	with tf.device('/gpu:0'):

		x_noise = tf.nn.dropout(x, keep_prob) #Apply noise to expression values
		for g in g_sets:
			n_feat                     = g_sizes[g]

			main_dict["W0_e_%s"%g] 	   = weight_variable([ n_feat ,1])

			main_dict["loc_m_%s"%g]    = tf.transpose(tf.concat([[x_noise[:,i]] for i in g_loc[g]], 0)) #[n_genes, batch]
			main_dict["h0_e_%s"%g]     = tf.matmul(main_dict["loc_m_%s"%g], main_dict["W0_e_%s"%g]) #no bias needed

		main_dict["h0_e"]  = dense_batch_lrelu(tf.concat([main_dict["h0_e_%s"%g] for g in g_sets], 1),
											   phase_train, scope="h0_e") #bn includes bias
		main_dict["h0d_e"] = tf.nn.dropout(main_dict["h0_e"], keep_prob) # Dropout

		# Additional layers after gene set layer
		l = -1 # In case there are no additional layers
		for l in xrange(len(layers)-1):
			main_dict["W%s_e"%(l+1)]   = weight_variable([ layers[l] , layers[l+1]])

			main_dict["h%s_e"%(l+1)]   = dense_batch_lrelu(tf.matmul(main_dict["h%sd_e"%l], main_dict["W%s_e"%(l+1)]), 
															phase_train, scope="h%s_e"%(l+1)) 
			main_dict["h%sd_e"%(l+1)]  = tf.nn.dropout(main_dict["h%s_e"%(l+1)], keep_prob)

	# DECODING
	with tf.device('/gpu:0'):

		main_dict["h0d_d"] = main_dict["h%sd_e"%(l+1)]
			
		# Decode hidden
		for k in xrange(len(rev_layers)-1):
			main_dict["W%s_d"%(k+1)]  = weight_variable([rev_layers[k], rev_layers[k+1]])

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

	sess        = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=False))
	sess.run(tf.global_variables_initializer())

	# Validation checks
	best_test      = np.Inf
	patience       = 0 # Epoch-based
	max_patience   = 50 # Epoch-based
	max_training   = 8000 # Maximum number of epochs to train for
	count_batch    = 0 # To get tabs of batch number
	count          = 0 # Number of epochs
	check          = 20 # Batch size during cost checking
	train_checks   = int(math.ceil(train_feat.shape[0] / float(check) )) #n_samples/check
	test_checks    = int(math.ceil(test_feat.shape[0] / float(check) )) #n_samples/check
	n_batches      = train_feat.shape[0] / batch_size # n_samples/batch_size

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

			test_error     = [cost.eval(feed_dict =
												{x:test_feat[i*check : (i+1)*check,],
												 keep_prob:1.0, phase_train: False})
							  for i in xrange(test_checks)]
		  	test_error     = np.mean(test_error)

		  	# Log errors
		  	log_table.append([count, batch_size, keepprob, train_error, test_error])

			# Update validation
			if test_error < best_test:
				best_test=test_error
				patience=0

				# UPDATE: Evaluate features for best model
				train_encoded = coded_feat.eval(feed_dict = {x:train_feat, keep_prob:1.0, phase_train:False})
				test_encoded  = coded_feat.eval(feed_dict = {x:test_feat, keep_prob:1.0, phase_train:False})
			else:
				patience+=1
		else:
			count_batch+=1

	# Clean up features
	train_encoded = pd.DataFrame(train_encoded,
								 columns=[str(i) for i in xrange(layers[-1])],
								 index=list(cell_train.index))
	test_encoded  = pd.DataFrame(test_encoded,
								 columns=[str(i) for i in xrange(layers[-1])],
								 index=list(cell_test.index))

	# Clean up log
	log_table     = pd.DataFrame(log_table, 
								 columns=["Epoch","bs","keep_prob", "train_error", "test_error"])
	sess.close()
	print("Done modeling GSEA-DAE")

	# Return
	return cell_train, cell_test, log_table

def gsea_autoencoder(cell_train, cell_test, in_folder, gsea="c2setcover", filter=int(10), arch="2"):
	# Builds autoencoder model for cell_train and applies to cell_test
	# gsea: c2, c2_cancer, c2setcover...

	# Load gsea data
	gsea = load_gsea_updated(gsea)

	# Build training feature out of cell training data
	gsea, cell_train, cell_test = gene_filter(gsea, cell_train, cell_test, filter=int(10))

	# Build autoencoded features
	cell_train, cell_test, log  = dae(cell_train, cell_test, gsea, arch)

	# Return
	return cell_train, cell_test, log