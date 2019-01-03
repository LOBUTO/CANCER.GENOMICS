# Set of machine learning models to predict pgms
import numpy as np
import pandas as pd
import random, math, gc, sys
from sklearn import preprocessing

# sys.path.insert(0, '/home/ubuntu/GIT')
sys.path.insert(0, '/tigress/zamalloa/GIT')
from paper_1 import parse_features_v2

def arch_layers(start, arch):
	arch   = arch.split("_")

	layers = [start]
	for i in arch:
		layers.append(start/int(i))

	return layers

def scale_standard_multiple(x, y, z=None, fix_index=False):
	# Will zero-center and divide by standard dev.
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per column across dataset standard: multiple")

	if fix_index==True:
		x_index = list(x.index)
		y_index = list(y.index)

	scaler = preprocessing.StandardScaler()
	x = scaler.fit_transform(x)

	if (y is not None):
		y = scaler.transform(y)

	if (z is not None):
		z_index = list(z.index)
		z = scaler.transform(z)

		if fix_index==True:
			x, y, z = pd.DataFrame(x, index=x_index), pd.DataFrame(y, index=y_index), pd.DataFrame(z, index=z_index)
		return x, y, z
	else:
		if fix_index==True:
			x, y = pd.DataFrame(x, index=x_index), pd.DataFrame(y, index=y_index)
		return x, y

def scale_0_1_multiple(x, y, z=None):
	# Will scale between 0-1 by column
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per gene across dataset 0-1: multiple")
	min_max_scaler = preprocessing.MinMaxScaler()
	x = min_max_scaler.fit_transform(x)

	if (y is not None):
		y = min_max_scaler.transform(y)

	if (z is not None):
		z = min_max_scaler.transform(z)
		return x, y, z
	else:
		return x, y

def process_test_target_features(cell_test, drug_test, target):
	# Processes target information for testing dataset
	
	target    = pd.read_csv(target, sep="\t", names=["Drug", "Sample"])

	# Consolidate
	samples   = list(set(target.Sample).intersection(list(cell_test.index)))

	# Do we need drug features?
	if drug_test is not None:
		drugs     = list(set(target.Drug).intersection(list(drug_test.index)))

		target    = (target
				    .query("Drug in @drugs")
				    .query("Sample in @samples"))

		# Build feature space
		data_feat = pd.concat([cell_test.loc[target.Sample,].reset_index(drop=True),
							   drug_test.loc[target.Drug,].reset_index(drop=True)],
							   axis=1)
	else:
		print("No drug features used!!!")
		target    = (target
					.query("Sample in @samples"))

		# Build feature space
		data_feat = cell_test.loc[target.Sample,].reset_index(drop=True)

	data_feat = data_feat.reset_index(drop=True)

	# Return
	print("Found {} drug-sample pairs to predict".format(data_feat.shape[0]))
	return target, data_feat

def keras_mlp(train_feat, train_target, test_feat, test_target,
			  batch_size, slr, keepprob, arch, l2_reg, patience=15, v_split=0.1, tuning=True, best_epoch=None):
	import random, math, gc
	from sklearn.metrics import mean_squared_log_error,mean_squared_error, r2_score,mean_absolute_error

	# Deep Learning Libraries
	import tensorflow as tf
	from keras.models import Sequential, load_model
	from keras.layers import Dense, Dropout, Flatten
	from keras.layers import Conv2D, MaxPooling2D, BatchNormalization
	from keras.layers.advanced_activations import PReLU
	from keras.optimizers import Adam,Nadam,SGD,Adagrad,Adadelta,RMSprop
	from keras.preprocessing.image import ImageDataGenerator
	from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
	from keras.utils import to_categorical
	from keras import backend as K
	from keras import metrics, regularizers, initializers
	# NOTE: Requires latest Keras version 2.2.4 (To extract best model at EarlyStopping)

	# Prep data
	train_feat, test_feat  = np.asarray(train_feat), np.asarray(test_feat)
	train_target = np.asarray(train_target.value)

	# Model as regression
	# Standardize initial inputs to 0-1 range
	# train_feat, test_feat = scale_0_1_multiple(train_feat, test_feat)
	# train_feat, test_feat = scale_standard_multiple(train_feat, test_feat)

	# Architecture
	init_features = train_feat.shape[1]
	print("init arch: ", init_features, arch)
	layers        = arch_layers(init_features, arch)
	print(layers)

	# Build KERAS Fully connected network
	SEED  = 1234
	model = Sequential()

	model.add(Dense(layers[1], input_dim=layers[0],
					kernel_regularizer=regularizers.l2(l2_reg), kernel_initializer = initializers.he_uniform(seed=SEED)))
	model.add(PReLU())
	model.add(BatchNormalization())
	model.add(Dropout(keepprob))

	for l in layers[2:]:
		model.add(Dense(l, kernel_regularizer=regularizers.l2(l2_reg), kernel_initializer = initializers.he_uniform(seed=SEED)))
		model.add(PReLU())
		model.add(BatchNormalization())
		model.add(Dropout(keepprob))

	# Add the output layer and compile
	model.add(Dense(1, kernel_initializer = initializers.he_uniform(seed=SEED)))

	# Optimizer
	optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)
	# optimizer = Adam(lr=slr, beta_1=0.9, beta_2=0.999)

	# Compile
	model.compile(optimizer=optimizer, loss="mse", metrics=["mse"])

	print(K.tensorflow_backend._get_available_gpus())
	sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))

	# Running model depending on need
	if tuning==True:
		early_stopping = EarlyStopping(monitor='val_mean_squared_error', patience=patience, restore_best_weights=True)

		print("Hyperparameter tunning....")
		print(model.summary())

		results     = model.fit(x=train_feat,
								y=train_target,
								batch_size=batch_size, epochs=4000, validation_split=v_split, verbose=2,
								shuffle=True,
								callbacks=[early_stopping])

		# Get stats
		best_metric = np.min(results.history["val_mean_squared_error"])
		train_loss  = results.history["mean_squared_error"][results.history["val_mean_squared_error"].index(best_metric)]
		best_epoch  = np.argmin(results.history["val_mean_squared_error"])
		log_table   = pd.DataFrame([[batch_size, keepprob, train_loss, best_metric]], 
									 columns=["bs","keep_prob", "train_error", "test_error"])

		print("Best Epoch found at {}".format(best_epoch))
		# Clean up
		K.clear_session()
		del model
		del results
		gc.collect()

		return best_epoch, log_table

	else: 
		print("Using Best Epoch found at {}".format(best_epoch))
		print("Building model....")
		print(model.summary())

		results     = model.fit(x=train_feat, y=train_target,
								batch_size=batch_size, epochs=best_epoch+1, verbose=2, shuffle=True)

		# Predict
		predictions = model.predict(test_feat)
		test_target = test_target.assign(prediction=predictions)
		
		# Clean up
		K.clear_session()
		del model
		del results
		gc.collect()

		return test_target

def drug_scale_nonbinary(drug_train, drug_test, fix_index=True):
	# Scales drug features, but leaves binary features, such as smiles features alone
	print("Scaling non-smiles drug features")

	all_cols      = list(drug_train.columns)
	smiles_labels = filter(lambda x: "s_" in x, all_cols)
	rest_labels   = filter(lambda x: "s_" not in x, all_cols)

	print("max before fitting: {}".format(np.max(np.max(drug_train[rest_labels]))))

	fixed_train, fixed_test = scale_standard_multiple(drug_train[rest_labels],
													  drug_test[rest_labels],
													  fix_index=fix_index)

	drug_train    = pd.concat([drug_train[smiles_labels], fixed_train], axis=1)
	drug_test     = pd.concat([drug_test[smiles_labels], fixed_test], axis=1)

	print("max after fitting: {}".format(np.max(np.max(drug_train))))
	return drug_train, drug_test

def cell_scale_nonbinary(cell_train, cell_test, fix_index=True):
	print("Scaling non-variant cell features")	

	all_cols      = list(cell_train.columns)
	var_labels    = filter(lambda x: "_var" in x, all_cols)
	rest_labels   = filter(lambda x: "_var" not in x, all_cols)

	print("max before fitting: {}".format(np.max(np.max(cell_train[rest_labels]))))

	fixed_train, fixed_test = scale_standard_multiple(cell_train[rest_labels],
													  cell_test[rest_labels],
													  fix_index=fix_index)

	cell_train    = pd.concat([cell_train[var_labels], fixed_train], axis=1)
	cell_test     = pd.concat([cell_test[var_labels], fixed_test], axis=1)

	print("max after fitting: {}".format(np.max(np.max(cell_train))))
	return cell_train, cell_test

def deep_pgms(target_table, cell_train, drug_train, cell_test, drug_test, target, 
			  arch="16_16_16_16", l2_reg=0.001, batch_size=200, slr=0.0001, 
			  keepprob=0.5, patience=15, cell_scaling=True, drug_scaling=True,
			  v_split=0.1):
	# Deep learning model to predict testing pgm modules

	# Do wee need to scale cell features
	if cell_scaling==True:
		cell_train, cell_test = cell_scale_nonbinary(cell_train, cell_test, fix_index=True)
	if drug_scaling==True:
		drug_train, drug_test = drug_scale_nonbinary(drug_train, drug_test, fix_index=True)

	# Process training features
	target_table, train_features = parse_features_v2(target_table, cell_train, drug_train)

	# Process test features
	test_table, test_features    = process_test_target_features(cell_test, drug_test, target)

	# import cPickle as pickle
	# with open("temp_feat.pickle", 'wb') as f:
	# 	pickle.dump([target_table, train_features, test_table, test_features], f)

	# Model tuning
	best_epoch, log_table = keras_mlp(train_features, target_table,
										test_features, test_table,
										batch_size, slr, keepprob, arch, l2_reg, patience, v_split, tuning=True)

	# Model building
	prediction   		  = keras_mlp(train_feat=train_features, train_target=target_table,
									  test_feat=test_features, test_target=test_table,
									  batch_size=batch_size, slr=slr, keepprob=keepprob, arch=arch, 
									  l2_reg=l2_reg, tuning=False, best_epoch=best_epoch)

	return prediction, log_table