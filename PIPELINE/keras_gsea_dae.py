import pandas as pd
import numpy as np
import tensorflow as tf
import random, math
from sklearn import preprocessing

# Set of Functions for Gene-set autoencoder
def load_gsea_updated(gsea_choice, in_folder):
	g_set = gsea_choice.split("_")
	
	all_g = []

	for g in g_set:
		file_in = "{}GSEA_FILES/{}_sets".format(in_folder, g)
		all_g.append(pd.read_csv(file_in, sep="\t"))

	all_g = pd.concat(all_g)

	return all_g

def gene_filter(gsea, train_exp, test_exp, filter):
	# Filters gsea and expression matrix ([samples x genes])

	# Find common genes
	common_genes   = set(train_exp.columns).intersection(list(test_exp.columns))
	common_genes   = list(set(gsea.genes).intersection(common_genes))

	# Apply threshold
	gsea           = gsea.loc[gsea.genes.isin(common_genes)]
	filtered_gsea  = (gsea.Gene_set.value_counts()
					  .reset_index(name="count")
					  .query("count<500")
					  .query("count>@filter")
					  .rename(columns={"index":"Gene_set"}))
	filtered_gsea  = gsea.loc[gsea.Gene_set.isin(list(filtered_gsea.Gene_set))]

	# Clean up
	filtered_genes = list(set(filtered_gsea.genes))
	train_exp      = train_exp[filtered_genes]
	test_exp       = test_exp[filtered_genes]

	# Return
	print("Total feature space: {}".format(train_exp.shape[1]))
	return filtered_gsea, train_exp, test_exp

def scale_standard_multiple(x, y, z=None):
	# Will zero-center and divide by standard dev.
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per gene across dataset standard: multiple")
	scaler = preprocessing.StandardScaler()
	x = scaler.fit_transform(x)

	if (y is not None):
		y = scaler.transform(y)

	if (z is not None):
		z = scaler.transform(z)
		return x, y, z
	else:
		return x, y

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

def arch_layers(start, arch):
	arch   = arch.split("_")

	layers = [start]
	for i in arch:
		layers.append(start/int(i))

	return layers

def dae(cell_train, cell_test, gsea, arch, keepprob, slr, batch_size, init_decomp, self_valid=True):
	from keras.layers import Input, Dense, Dropout, BatchNormalization, concatenate, Lambda
	from keras.layers.advanced_activations import PReLU
	from sklearn.metrics import mean_squared_log_error,mean_squared_error, r2_score,mean_absolute_error
	from keras.optimizers import Nadam, Adam,SGD,Adagrad,Adadelta,RMSprop
	from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
	from keras.models import Model
	from tensorflow import gather
	from keras import backend as K
	from keras import metrics, regularizers, initializers
	import gc
	# Trains autoencoder using tensorflow 

	# Prepare data
	g_sets        = list(set(gsea.Gene_set))
	all_genes     = list(set(gsea.genes))
	exp_cols      = list(cell_train.columns)
	index_list    = []
	for i in g_sets:
		genes = list(gsea.query("Gene_set==@i").genes)
		index = [exp_cols.index(i) for i in genes]
		index_list.append(index)

	train_feat    = np.asarray(cell_train) #[samples, genes]
	test_feat     = np.asarray(cell_test)
	log_table     = []

	# Standardize initial inputs to 0-1 range
	# train_feat, test_feat = scale_0_1_multiple(train_feat, test_feat)
	train_feat, test_feat = scale_standard_multiple(train_feat, test_feat)
	print("Number of gene sets used for parallel encoding: {}".format(len(g_sets)))

	# Specify architecture
	layers        = arch_layers(len(index_list)*init_decomp, arch)
	rev_layers    = [i for i in reversed(layers)]
	print("layers: ", layers)
	print("rev_layers:", rev_layers)

	# BUILD KERAS MODEL - API
	SEED 		= 1234
	recon_feat  = len(all_genes)
	input_shape = Input(shape=(train_feat.shape[1],))

	# Split data based on index list
	for il in xrange(len(index_list)):
		vars()["x_%s"%il] = Lambda(lambda x: gather(x, index_list[il], axis=1 ))(input_shape)

	# Learn each individually
	for il in xrange(len(index_list)):
		vars()["x_%s_prel"%il] = Dense(init_decomp, activation=PReLU(), 
										kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["x_%s"%il])

	# Merge for first layer
	encoded_merge                = concatenate([vars()["x_%s_prel"%il] for il in xrange(len(index_list))])
	encoded_bn                   = BatchNormalization()(encoded_merge)
	vars()["encode_%s_drop"%(0)] = Dropout(keepprob)(encoded_bn)

	# Add additional encoded layers [200,100]
	for l in xrange(len(layers)-1):
		vars()["encode_%s_prel"%(l+1)] = Dense(layers[(l+1)], activation=PReLU(),
											kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["encode_%s_drop"%l])
		vars()["encode_%s_bn"%(l+1)]   = BatchNormalization()(vars()["encode_%s_prel"%(l+1)])
		vars()["encode_%s_drop"%(l+1)] = Dropout(keepprob)(vars()["encode_%s_bn"%(l+1)])

	# DECODING [100,200]
	# Decoding first layer
	vars()["decode_0_prel"] = Dense(rev_layers[1], activation=PReLU(),
									kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["encode_%s_drop"%(l+1)])
	vars()["decode_0_bn"]   = BatchNormalization()(vars()["decode_0_prel"])
	vars()["decode_0_drop"] = Dropout(keepprob)(vars()["decode_0_bn"])
	
	# Decoding additional layers
	d=-1
	for d in xrange(len(rev_layers)-2):
		vars()["decode_%s_prel"%(d+1)] = Dense(rev_layers[(d+2)], activation=PReLU(),
											kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["decode_%s_drop"%d])
		vars()["decode_%s_bn"%(d+1)]   = BatchNormalization()(vars()["decode_%s_prel"%(d+1)])
		vars()["decode_%s_drop"%(d+1)] = Dropout(keepprob)(vars()["decode_%s_bn"%(d+1)])

	# Reconstruct
	vars()["last"] = Dense(recon_feat, activation="sigmoid")(vars()["decode_%s_drop"%(d+1)])

	# Map the autoencoder
	autoencoder = Model(input_shape, vars()["last"])

	# Create encoding layer (Latent space)
	encoder     = Model(input_shape, vars()["encode_%s_prel"%(l+1)])

	# Optimizer
	optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)

	# Compile
	autoencoder.compile(optimizer=optimizer, loss="mean_squared_error", metrics=["mse"])
	
	# Run model
	early_stopping = EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True)
	
	print(autoencoder.summary())

	# How are we approaching validation - On a split of train set or validating on external set
	if self_valid==True:
		results        = autoencoder.fit(x=train_feat, y=train_feat, validation_split=0.4, #validation_data=(test_feat, test_feat),
									batch_size=batch_size, verbose=2, epochs=500,
									callbacks=[early_stopping])
	else:
		results        = autoencoder.fit(x=train_feat, y=train_feat, validation_data=(test_feat, test_feat),
									batch_size=batch_size, verbose=2, epochs=500,
									callbacks=[early_stopping])
	
	# Obtain prediction
	train_output   = encoder.predict(train_feat)
	test_output    = encoder.predict(test_feat)

	# Clean up
	K.clear_session()
	del autoencoder
	gc.collect()
	print("Done modeling GSEA-DAE")

	# Return
	return train_output, test_output, results

def gsea_autoencoder(cell_train, cell_test, in_folder, gsea="c2setcover", filter_th=int(10), arch="2", labeller="_gs",
					init_decomp=5, self_valid=True):
	# Builds autoencoder model for cell_train and applies to cell_test
	# gsea: c2, c2_cancer, c2setcover...
	print("GSEA Autoencoded Features...")

	# Produce hyperparameters
	keepprob        = float(0.5)
	slr             = 0.001
	batch_size      = 20
	init_decomp     = init_decomp # To how many neurons is each gene_set originally being decomposed to

	# Load gsea data
	gsea = load_gsea_updated(gsea, in_folder)

	# Build training feature out of cell training data
	gsea, cell_train, cell_test = gene_filter(gsea, cell_train, cell_test, filter=filter_th)

	# Obtain labels
	train_labels  = list(cell_train.index)
	test_labels   = list(cell_test.index)

	# Build autoencoded features
	cell_train, cell_test, results  = dae(cell_train, cell_test, gsea, arch, keepprob, slr, 
										  batch_size, init_decomp, self_valid=self_valid)

	# Get stats
	best_metric   = np.min(results.history["val_loss"])
	train_loss    = results.history["loss"][results.history["val_loss"].index(best_metric)]
	log_table     = pd.DataFrame([[batch_size, keepprob, train_loss, best_metric]], 
								 columns=["bs","keep_prob", "train_error", "test_error"])

	# Clean up and return
	cell_train    = pd.DataFrame(cell_train, index=train_labels,
								columns = [str(i) + labeller for i in xrange(cell_train.shape[1])])
	cell_test     = pd.DataFrame(cell_test, index=test_labels,
								columns = [str(i) + labeller for i in xrange(cell_test.shape[1])])

	# Manual loading
	# cell_train = pd.read_csv("/home/ubuntu/DAE_FILES/gdsc_ctrp_train_gsea_c2setcovertanimoto_dae.txt", sep="\t", index_col=0)
	# cell_test  = pd.read_csv("/home/ubuntu/DAE_FILES/gdsc_ctrp_test_gsea_c2setcovertanimoto_dae.txt", sep="\t", index_col=0)
	# log_table  = pd.read_csv("/home/ubuntu/DAE_FILES/gdsc_ctrp_log_c2setcovertanimoto_dae.txt", sep="\t")

	# Return
	return cell_train, cell_test, log_table