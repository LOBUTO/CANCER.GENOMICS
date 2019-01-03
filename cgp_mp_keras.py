import numpy as np
import pandas as pd
import math 
import gc
import itertools
import random
import sys
## Libraries

from sklearn.metrics import mean_squared_log_error,mean_squared_error, r2_score,mean_absolute_error
from imblearn.over_sampling import SMOTE, ADASYN
import tensorflow as tf

# Deep Learning Libraries
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv2D, MaxPooling2D, BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.optimizers import Adam,Nadam,SGD,Adagrad,Adadelta,RMSprop
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
from keras.utils import to_categorical
from keras import backend as K
from keras import metrics

# Import functions
sys.path.insert(0, '/tigress/zamalloa/GIT')

from paper_1 import *

def keras_mlp(train_feat , valid_feat , train_target , valid_target, 
              arch, slr, keepprob, problem, pred_file, batch_size, epochs):

    # Params
    batch_size = batch_size
    epochs     = epochs

    # Standardize initial inputs to 0-1 range
    train_feat, valid_feat = scale_0_1_multiple(train_feat, valid_feat)
    print("SHAPES!!!", train_feat.shape, valid_feat.shape)

    if problem=="classification":
        print("Classification problem")
        
        # Get class weights
        ones    = np.count_nonzero(train_target)
        zeros   = train_target.shape[0] - ones
        zero_weight = float(zeros)/ones
        
        # class_weight = {0:1.0, 1:zero_weight}
#         class_weight = {0:zero_weight, 1:1.0}
        class_weight = {0:1.0, 1:2.5}
        print(class_weight)
        
        # One-hot encode labels (2-Classes)
        y_train = to_categorical(train_target, num_classes=2)
        y_valid = to_categorical(valid_target, num_classes=2)

    elif problem=="regression":
        y_train      = train_target
        y_valid      = valid_target
        class_weight = None

    # Architecture
    init_features = train_feat.shape[1]
    print("init arch: ", init_features, arch)
    layers        = arch_layers(init_features, arch)
    print(layers)

    # Build KERAS Fully connected network
    model = Sequential()
    
    
    model.add(Dense(layers[1], activation=PReLU(), input_dim=layers[0]))
    model.add(BatchNormalization())
    model.add(Dropout(keepprob))    

    for l in layers[2:]:
        model.add(Dense(l, activation=PReLU()))
        model.add(BatchNormalization())
        model.add(Dropout(keepprob))

    # Add the output layer and compile
    if problem=="classification": 
        model.add(Dense(2, activation='softmax'))
        
        # Optimizer
        optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)

        # Compile
        model.compile(optimizer=optimizer, loss="binary_crossentropy", 
                      metrics=[metrics.binary_crossentropy])

    elif problem=="regression":
        model.add(Dense(1))

        # Optimizer
        optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)

        # Compile
        model.compile(optimizer=optimizer, loss="mse", metrics=["mse"])
    
    print(model.summary())

    # Run model
#     reduce_lr      = LearningRateScheduler(lambda x: 1e-3 * 0.9 ** x)
    early_stopping = EarlyStopping(monitor='val_loss', patience=20)
#     results   = model.fit(train_feat, y_train, batch_size=batch_size, epochs=epochs,
#                           validation_data=(valid_feat, y_valid), verbose=2,
#                           callbacks=[reduce_lr], class_weight=class_weight)

#     model.fit(x=np.concatenate((train_feat, valid_feat), axis=0), 
#               y=np.concatenate((y_train, y_valid), axis=0), 
#               batch_size=batch_size, epochs=5, validation_split=0.1, verbose=2,
#               callbacks=[reduce_lr], class_weight=class_weight)
    results   = model.fit(x=train_feat, 
                          y=y_train, 
                          batch_size=batch_size, epochs=200, validation_split=0.1, verbose=2,
                          callbacks=[early_stopping], class_weight=class_weight)
    
    # Clean up
#     K.clear_session()
    
    return results, model, early_stopping

def run(exp_arch, exp_gfilter, drug_arch, keepprob, arch, slr, target_features, exp_target,
        encoded_target, encoded_target_arch,
        string_th, gdsc_variants, giant_features, problem, batch_size, epochs, drug):
    
    home_folder         = "/tigress/zamalloa/"
    out_folder          = "{}CGP_FILES/RESULTS/".format(home_folder)
     
    # Load data
    lobico, exp_data, mut_data, drug_data, lobico_ori = load_data(exp_arch, exp_gfilter, drug_arch, problem, in_folder=home_folder)

    # Build feature space
    lobico, data_feat   = parse_features(lobico, exp_data, mut_data, drug_data, target_features, exp_target, 
                                        encoded_target, encoded_target_arch, string_th, gdsc_variants, giant_features, 
                                        in_folder=home_folder)

    # Prep output files (for all folds)
    pred_file             = prep_output_files_keras(out_folder, "drugout", exp_arch, exp_gfilter, drug_arch, 
                                            keepprob, arch, slr, 
                                            target_features, exp_target, encoded_target, encoded_target_arch, string_th, 
                                            gdsc_variants, giant_features, problem, drug)
    
    # Split into training and validation
    x_train , x_valid , y_train , y_valid, true_values = drug_feat_target_split(lobico, data_feat, drug)

    # Run model per cell
    results, model, early_stopping = keras_mlp(x_train , x_valid , y_train , y_valid,
                                               arch, slr, keepprob, problem,
                                               pred_file, batch_size, epochs)
    best_metric = np.min(results.history["val_loss"]) #Given that accuracy is the early stopping metric

    # Predict
    cells       = lobico.loc[lobico.Compound==drug]["cell_name"].values
    print("Predicting...")
    probabilities   = model.predict(x_valid)

    if problem=="classification":
        with open(pred_file, "a") as f:
            for i in xrange(len(true_values)):
                f.write("\n%s\t%s\t%s\t%s\t%s\t%s"%(drug, cells[i], true_values[i], probabilities[i,0], probabilities[i,1], best_metric))

    elif problem=="regression":
        with open(pred_file, "a") as f:
            for i in xrange(len(true_values)):
                f.write("\n%s\t%s\t%s\t%s\t%s"%(drug, cells[i], true_values[i], probabilities[i][0], best_metric))
        
        # Clean up
        del model
        del results
        gc.collect()

def run_wrapper():
    # Arguments
    drug                = sys.argv[1]
    exp_arch            = "2" # The expression features dae (cgp_encode.py) results
    exp_gfilter         = int(10) # The expression g_filter used (cgp_encode.py)
    drug_arch           = "8_16" # "pubchem_arch_8_16" #"8_16" 
    arch                = "1_4"
    keepprob            = float(0.5)
    slr                 = 0.00001 #0.0000001
    target_features     = True
    exp_target          = False #If target features true, then we will use regression features, otherwise binary features
    encoded_target      = True
    encoded_target_arch = "2_16" #Only used if encoded_target==True
    string_th           = 550 #Irrelevant if both target_features==False & exp_target==False
    gdsc_variants       = True
    giant_features      = True
    problem             = "classification" #regression
    batch_size          = 200
    epochs              = 20
    
    drug                = Function_drug_name_zhang_to_gdsc(drug, "ztg")
    run(exp_arch, exp_gfilter, drug_arch, keepprob,
        arch, slr, target_features, exp_target, encoded_target, encoded_target_arch,
        string_th, gdsc_variants, giant_features, problem, batch_size, epochs,
        drug)
    print("DONE")

if __name__ == "__main__":
    run_wrapper()