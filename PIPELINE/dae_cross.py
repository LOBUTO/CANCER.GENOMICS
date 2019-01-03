import numpy as np
import pandas as pd
import random, math, gc
from sklearn import preprocessing

def arch_layers(start, arch):
	arch   = arch.split("_")

	layers = [start]
	for i in arch:
		layers.append(start/int(i))

	return layers

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

def cross_autoencoder(train, train_noisy, test, encode, decode, keepprob, slr, batch_size, tuning=True, best_epoch=None):
	from keras.layers import Input, Dense, Dropout, BatchNormalization
	from keras.layers.advanced_activations import PReLU
	from sklearn.metrics import mean_squared_log_error,mean_squared_error, r2_score,mean_absolute_error
	from keras.optimizers import Nadam, Adam,SGD,Adagrad,Adadelta,RMSprop
	from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
	from keras.models import Model
	from keras import backend as K
	from keras import metrics, regularizers, initializers
	import gc
	# NOTE: Requires latest Keras version 2.2.4 (To extract best model at EarlyStopping)

	print("Encoding layers")
	
	# Start
	SEED        = 1234
	input_shape = Input(shape=(train.shape[1],))
	print(encode)
	e           = encode[0] # In case single layer
	vars()["encoded_%s"%e] = Dense(encode[0], activation=PReLU(), kernel_initializer = initializers.he_uniform(seed=SEED))(input_shape)
	vars()["e_batch_%s"%e] = BatchNormalization()(vars()["encoded_%s"%e])
	vars()["e_drop_%s"%e]  = Dropout(keepprob)(vars()["e_batch_%s"%e])

	# Encoding layers
	prev_layer = encode[0]
	print("prev_layer: ", prev_layer)
	for e in encode[1:]:
		print(e)
		vars()["encoded_%s"%e] = Dense(e, activation=PReLU(), kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["e_drop_%s"%prev_layer])
		vars()["e_batch_%s"%e] = BatchNormalization()(vars()["encoded_%s"%e])
		vars()["e_drop_%s"%e]  = Dropout(keepprob)(vars()["e_batch_%s"%e])
		prev_layer = e

	# Decoding layers
	vars()["d_drop_%s"%prev_layer] = vars()["e_drop_%s"%e]
	d = prev_layer  # In case single layer
	for d in decode[:-1]:
		print(d, vars()["d_drop_%s"%prev_layer])
		vars()["decoded_{}".format(d)]  = Dense(d, activation=PReLU(), kernel_initializer = initializers.he_uniform(seed=SEED))(vars()["d_drop_{}".format(prev_layer)])
		vars()["d_batch_%s"%d] = BatchNormalization()(vars()["decoded_%s"%d])
		vars()["d_drop_%s"%d]  = Dropout(keepprob)(vars()["d_batch_%s"%d])
		prev_layer = d

	vars()["last"] = Dense(decode[-1], activation="sigmoid")(vars()["d_drop_%s"%d])

	# Map the autoencoder
	autoencoder = Model(input_shape, vars()["last"])

	# Create encoding layer
	encoder     = Model(input_shape, vars()["encoded_%s"%e])

	# Optimizer
	optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)

	# Compile
	autoencoder.compile(optimizer=optimizer, loss="mean_squared_error", metrics=["mse"])

	# Run model depending on need
	if tuning==True:
		early_stopping = EarlyStopping(monitor='val_loss', patience=100, restore_best_weights=True)

		print("Hyperparameter tunning....")
		print(autoencoder.summary())
		results        = autoencoder.fit(x=train_noisy, y=train, validation_split=0.2, #validation_data=(test, test), 
									batch_size=batch_size, epochs=2000, verbose=2,
									callbacks=[early_stopping])
		best_metric    = np.min(results.history["val_loss"])
		train_loss     = results.history["loss"][results.history["val_loss"].index(best_metric)]
		best_epoch     = np.argmin(results.history["val_loss"])
		K.clear_session()

		return best_epoch, best_metric, train_loss

	else: 
		print("Best Epoch found at {}".format(best_epoch))
		print("Building model....")
		results        = autoencoder.fit(x=train_noisy, y=train, epochs=best_epoch+1,
									batch_size=batch_size, verbose=2)

		# Obtain prediction
		train_output   = encoder.predict(train)
		test_output    = encoder.predict(test)
		K.clear_session()

		return train_output, test_output

def cross_feature_dae(train, test, arch, labeller=""):
	# Builds a DAE model based on train, and regulates it based on test
	# Produces Autoencoded matrices for both inputs

	# Produce hyperparameters
	keepprob        = float(0.5)
	slr             = 0.0001
	batch_size      = 20
	noise           = 0.0

	# Produce architecture
	init_features = train.shape[1]
	layers        = arch_layers(init_features, arch) #[i, i/2, i/4]
	rev_layers    = [i for i in reversed(layers)] #[i/4, i/2, i]

	# Obtain labels
	train_labels  = list(train.index)
	test_labels   = list(test.index)

	# Standardize initial inputs to 0-1 range
	train, test   = scale_0_1_multiple(train, test) # NOTE: Since values have a min of Zero

	# Noise data
	train_noisy   = train + noise * np.random.normal(loc=0.0, scale=1.0, size=train.shape)
	train_noisy   = np.clip(train_noisy, 0., 1.) #Clip to range to 0-1

	# Model tuning
	best_epoch, best_metric, train_loss   = cross_autoencoder(train, train_noisy, test,
												layers[1:], rev_layers[1:], keepprob, slr, batch_size,
												tuning=True)
	# Model building
	train, test    						  = cross_autoencoder(train, train_noisy, test,
												layers[1:], rev_layers[1:], keepprob, slr, batch_size,
												tuning=False, best_epoch=best_epoch)

	# Get stats table
	log_table     = pd.DataFrame([[batch_size, keepprob, train_loss, best_metric]], 
								 columns=["bs","keep_prob", "train_error", "test_error"])

	# Clean up and return
	train         = pd.DataFrame(train, index=train_labels,
								columns = [str(i) + labeller for i in xrange(train.shape[1])])
	test          = pd.DataFrame(test, index=test_labels,
								columns = [str(i) + labeller for i in xrange(test.shape[1])])

	return train, test, log_table