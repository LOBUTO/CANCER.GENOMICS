# MLP to predict all cells at once for each drug compound using all other information

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
import itertools

# Import functions
sys.path.insert(0, '/tigress/zamalloa/GIT')

from paper_1 import *
from tf_add import *

# Execute
out_folder    = "/tigress/zamalloa/CGP_FILES/RESULTS/"
exp_arch      = sys.argv[1] # "2" # The expression features dae (cgp_encode.py) results
exp_gfilter   = int(sys.argv[2]) # The expression g_filter used (cgp_encode.py)
drug_arch     = sys.argv[3] #"2" # The drug features dae (drug_dae.py) results
keepprob      = float(sys.argv[4])
arch          = sys.argv[5] # The architecture to use
drug          = sys.argv[6]
slr           = float(sys.argv[7]) #0.001 #0.0000001
# epsilon       = float(sys.argv[8])
target_features = bool(sys.argv[8]=="True")
exp_target      = bool(sys.argv[9]=="True") #If target features true, then we will use regression features, otherwise binary features

# Prep features and target data
exp_data      = exp_encode_parse("c2setcover", exp_gfilter, exp_arch, 0.5)
drug_data     = drug_encode_parse(drug_arch, 0.5)
# lobico        = process_paper_lobico("/tigress/zamalloa/PAPERS/LORIO/lobico.csv")

lobico        = cgp_act_post_process(pd.read_csv("/tigress/zamalloa/CGP_FILES/v17.3_fitted_dose_response.csv"),
                                    zscoring=False)
lobico, exp_data, drug_data = drug_exp_lobico_parse_v2(lobico, exp_data, drug_data, "ztg")

# Do we need to include target features
if target_features==True:
    print("Using drug target features")
    pp_hugo             = string_pp_to_hugo("/tigress/zamalloa/STRING/9606.protein.links.v10.5.txt",
                                            "/tigress/zamalloa/STRING/9606.protein.aliases.v10.5.txt")
    drug_target         = process_drug_target("/tigress/zamalloa/CGP_FILES/Screened_Compounds.csv")
    
    if exp_target==True:
        print("Using target expression features")
        drug_target         = string_target_binary_features(drug_target, pp_hugo, th=900, output="table")
        lobico, target_feat = string_target_expression_features(drug_target, "CGP_FILES/070818_cgp_exp.txt", lobico)
        
        data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                    drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                    target_feat.reset_index(drop=True)], axis=1)
    else:
        print("Using target binary features")
        drug_target         = string_target_binary_features(drug_target, pp_hugo, th=900, output="pivot")
        lobico, target_feat = target_lobico_parse(lobico, drug_target)

        data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                    drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                    target_feat.loc[lobico.Compound].reset_index(drop=True)], axis=1)
else:
    print("No target features used")
    data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                drug_data.loc[lobico.Compound,].reset_index(drop=True)], axis=1)

print(data_feat.head())
print("All data:", data_feat.shape)

# We are going to predict all cells at once per drug target
lobico        = lobico.reset_index(drop=True)

# Prep output files
log_file      = "{}cgpmlp_exp_drug_reg_{}_{}_{}_{}_{}_{}_{}_{}_{}_log".format(out_folder, drug, exp_arch, exp_gfilter, 
                                                                 drug_arch, keepprob, arch, slr, target_features, exp_target)
pred_file     = "{}cgpmlp_exp_drug_reg_{}_{}_{}_{}_{}_{}_{}_{}_{}_txt".format(out_folder, drug, exp_arch, exp_gfilter, 
                                                                 drug_arch, keepprob, arch, slr, target_features, exp_target)

with open(log_file, "w") as f:
    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%("Epoch", "bs", "keep_prob", "drug_name", "train_error",
                                              "valid_error", "test_error", "lr"))

with open(pred_file, "w") as f:
    f.write("Compound\tEpoch\tcell_name\tprediction\ttrue_value")

tf.reset_default_graph()

target_index = list(lobico.loc[lobico.Compound==drug].index)
cells        = list(lobico.loc[lobico.Compound==drug].cell_name)
print(len(target_index))

test_feat    = np.asarray(data_feat.loc[target_index,:])

#################### SPLIT BY RANDOM COMPOUND SET ####################
train_index, valid_index = data_split_drug_self(lobico, drug)
######################################################################

train_feat   = np.asarray(data_feat.loc[train_index])
valid_feat   = np.asarray(data_feat.loc[valid_index])

train_target = np.asarray(lobico.loc[train_index].value)
valid_target = np.asarray(lobico.loc[valid_index].value)
test_target  = np.asarray(lobico.loc[target_index].value)

# Standardize initial inputs to 0-1 range
train_feat, valid_feat, test_feat = scale_0_1_multiple(train_feat, valid_feat, test_feat)

print("SHAPES!!!", train_feat.shape, valid_feat.shape, test_feat.shape)

# Build graphs
init_features = train_feat.shape[1]
print("init arch: ", init_features, arch)
layers        = arch_layers(init_features, arch)
print(layers)

main_dict     = {}

# Placeholders
with tf.device('/gpu:0'):
    x           = tf.placeholder(tf.float32, shape=[None, init_features], name="x")
    y_          = tf.placeholder(tf.float32, shape=[None], name="y_")
    phase_train = tf.placeholder(tf.bool, name='phase_train')
    keep_prob   = tf.placeholder(tf.float32, name="keep_prob")
    pw          = tf.placeholder(tf.float32, name="pw")
    learning_rate = tf.placeholder(tf.float32, shape=[])

# Hidden Layers
with tf.device('/gpu:0'):

#             main_dict["h0d"] = tf.nn.dropout(x, 0.9) #Apply noise to features
    main_dict["h0d"] = x

    for l in xrange(len(layers)-1):
        main_dict["W%s"%(l+1)]  = weight_variable([ layers[l] , layers[l+1]])

        main_dict["h%s"%(l+1)]  = dense_lrelu_batch(tf.matmul(main_dict["h%sd"%l], main_dict["W%s"%(l+1)]), 
                                                    phase_train, scope="h%s%s"%(l+1, drug))
        main_dict["h%sd"%(l+1)] = tf.nn.dropout(main_dict["h%s"%(l+1)], keep_prob)

# Output layer
W_fcn  = weight_variable([layers[-1], 1])
b_fcn  = bias_variable([1])
y_fc   = tf.add(tf.matmul(main_dict["h%sd"%(l+1)], W_fcn),  b_fcn, name="y_fc")

# Calculate error with original features (denoised)
with tf.name_scope('cost'):
    cost  = tf.reduce_mean(tf.square(tf.transpose(y_fc)[0] - y_))
    
    prediction  = tf.transpose(y_fc)[0]
    true_target = y_

#################################################### RUNNING ########################################################
update_ops  = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
with tf.control_dependencies(update_ops):
    # Ensures that we execute the update_ops before performing the train_step
    # train_step    = tf.contrib.opt.NadamOptimizer(slr).minimize(cost)
    # train_step    = NadamOptimizer(slr).minimize(cost)
    train_step            = tf.train.AdamOptimizer(slr).minimize(cost)

sess          = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
sess.run(tf.global_variables_initializer())

# Parameters
batch_size    = 200

# Validation checks
best_valid     = np.Inf
patience       = 0 # Epoch-based
max_patience   = 100 # Epoch-based
max_training   = 500 # Maximum number of epochs to train for
count_batch    = 0 # To get tabs of batch number
count          = 0 # Number of epochs
check          = 200 # Batch size during cost checking
train_checks   = int(math.ceil(train_feat.shape[0] / float(check) )) #n_samples/check
valid_checks   = int(math.ceil(valid_feat.shape[0] / float(check) ))
n_batches      = train_feat.shape[0] / batch_size # n_samples/batch_size


print("n_batches: %s"%n_batches)
print("train_checks: %s"%train_checks)
while (patience < max_patience) and (count < max_training):
#         while (slr > 1e-10) and (count < max_training):
    
    batch_all    = random.sample(xrange(train_feat.shape[0]), batch_size)
    feat_batch   = train_feat[batch_all,]
    target_batch = train_target[batch_all]
    
    train_step.run(feed_dict={x: feat_batch, y_: target_batch,
                              keep_prob: keepprob, phase_train: True})#,
                             #learning_rate:slr}) #, pw:pos_weights})

    if count_batch == (n_batches-1): # Batch limit count, that is a whole epoch has passed
        count_batch = 0
        count+=1 # Epoch number increase

        train          = [sess.run([cost],
                                   feed_dict =
                                        {x:train_feat[i*check : (i+1)*check,],
                                         y_:train_target[i*check : (i+1)*check],
                                         keep_prob:1.0, phase_train: False})
                          for i in xrange(train_checks)]

        train_cost     = np.mean([i[0] for i in train])
        
        validation     = [sess.run([cost],
                                    feed_dict =
                                        {x:valid_feat[i*check : (i+1)*check,],
                                         y_:valid_target[i*check : (i+1)*check,],
                                         keep_prob:1.0, phase_train:False})
                          for i in xrange(valid_checks)]
        valid_cost     = np.mean([i[0] for i in validation])
        
        test_cost, pred, true = sess.run([cost, prediction, true_target],
                                         feed_dict={x:test_feat, y_:test_target, keep_prob:1.0, phase_train:False})
        # Write to log
        with open(log_file, "a") as f:
            f.write("\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(count, batch_size, keepprob, drug, 
                                                        train_cost, valid_cost, test_cost, slr))
        # print("count, train, valid, test, valid_roc: ", count, train_cost, valid_cost, test_cost, patience, slr)
        
        # Write predictions
        with open(pred_file, "a") as f:
            for m in xrange(len(cells)):
                f.write("\n%s\t%s\t%s\t%s\t%s"%(drug, count, cells[m], pred[m], true[m]))
        
        # Update validation
        if (valid_cost < best_valid) and (count>5):
            best_valid=valid_cost
            patience=0
            
        else:
            patience+=1
    else:
        count_batch+=1
sess.close()
print("Done modeling pack-cell MLP")

print("DONE")