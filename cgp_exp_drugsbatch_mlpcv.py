#### drug_dae.py
# MLP to predict all cells at once for each drug compound using all other information
# This is for cross-validation purposes, that is we hold a fifth of the drugs out, we train on the rest
# and we calculate and ROC per drug, we stop training after x-patience

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
from sklearn import metrics
from sklearn.metrics import roc_auc_score

# Import functions
sys.path.insert(0, '/tigress/zamalloa/GIT')

from paper_1 import *

def exp_encode_parse(gsea, g_filter, arch, dropkeep):
    in_folder = "/tigress/zamalloa/GSEA_FILES/RESULTS/"
    
    main_list = []
    for i in ["train", "valid", "test"]:
        in_file = "{}CGP_{}_{}_deepautoencoder_{}_{}_{}.txt".format(in_folder,gsea,g_filter,arch,dropkeep,i)
        
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

def drug_encode_parse(arch, dropkeep):
    in_folder = "/tigress/zamalloa/CGP_FILES/RESULTS/"
    
    main_list = []
    for i in ["train", "valid", "test"]:

        in_file = "{}{}_deepautoencoder_{}_{}_{}.txt".format(in_folder, "cgp_drugs", arch, dropkeep,i)
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

def load_data(exp_arch, exp_gfilter, drug_arch):
    # Data for ROC
    lobico_ori = process_paper_lobico("/tigress/zamalloa/PAPERS/LORIO/lobico.csv")
    lobico_ori["Compound"] = [Function_drug_name_zhang_to_gdsc(i, "ztg") for i in lobico_ori["Compound"]]

    # Prep data
    exp_data      = exp_encode_parse("c2setcover", exp_gfilter, exp_arch, 0.5)
    if drug_arch=="original":
        drug_data     = pd.read_csv("/tigress/zamalloa/PAPERS/ZHANG_2018/drug_feats_proc.txt", sep="\t", index_col=0)
    elif drug_arch=="pubchem_smiles":
        drug_data     = pd.read_csv("/tigress/zamalloa/CGP_FILES/pubchem_smiles.txt", index_col=0)
    else:
        drug_data     = drug_encode_parse(drug_arch, 0.5)
        
    lobico        = cgp_act_post_process(pd.read_csv("/tigress/zamalloa/CGP_FILES/v17.3_fitted_dose_response.csv"),
                                        zscoring=False)
    lobico, exp_data, drug_data = drug_exp_lobico_parse_v2(lobico, exp_data, drug_data, "ztg")

    return lobico, exp_data, drug_data, lobico_ori

def parse_features(lobico, exp_data, drug_data, target_features, exp_target, string_th):
    # Do we need to add target features??
    if target_features==True:
        print("Using drug target features")
        pp_hugo             = string_pp_to_hugo("/tigress/zamalloa/STRING/9606.protein.links.v10.5.txt",
                                                "STRING/9606.protein.aliases.v10.5.txt")
        drug_target         = process_drug_target("/tigress/zamalloa/CGP_FILES/Screened_Compounds.csv")

        if exp_target==True:
            print("Using target expression features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="table") #[Compound, target, value==1]
            lobico, target_feat = string_target_expression_features(drug_target, "CGP_FILES/070818_cgp_exp.txt", lobico)

            data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                        drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                        target_feat.reset_index(drop=True)], axis=1)
        else:
            print("Using target binary features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="pivot")
            lobico, target_feat = target_lobico_parse(lobico, drug_target)

            data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                        drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                        target_feat.loc[lobico.Compound].reset_index(drop=True)], axis=1)
        
    else:
        print("No target features used")
        data_feat      = pd.concat([exp_data.loc[lobico.cell_name,].reset_index(drop=True),
                                    drug_data.loc[lobico.Compound,].reset_index(drop=True)], axis=1)
    
    lobico        = lobico.reset_index(drop=True)
    print(data_feat.head())
    print("All data:", data_feat.shape)

    return(lobico, data_feat)

def prep_output_files(out_folder, exp_arch, exp_gfilter, drug_arch, keepprob, arch, slr, target_features, exp_target, string_th):
    # Prep output files
    log_file      = "{}cgpmlp_exp_drug_reg_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_log".format(out_folder, "drugscv", exp_arch, exp_gfilter, 
                                                                     drug_arch, keepprob, arch, slr, target_features, exp_target, string_th)
    pred_file     = "{}cgpmlp_exp_drug_reg_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_txt".format(out_folder, "drugscv", exp_arch, exp_gfilter, 
                                                                     drug_arch, keepprob, arch, slr, target_features, exp_target, string_th)

    with open(log_file, "w") as f:
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s"%("CFOLD", "Epoch", "bs", "keep_prob", "train_error",
                                                  "valid_error", "lr"))

    with open(pred_file, "w") as f:
        f.write("CFOLD\tEpoch\tCompound\tcell_name\tprediction\ttrue_value\n")

    return log_file, pred_file

def chunk_split(lobico, folds=5):
    # Split into n drug folds
    all_drugs   = list(set(lobico.Compound))
    drug_chunks = [list(i) for i in np.array_split(all_drugs, folds)]

    return all_drugs, drug_chunks

def mlp(lobico, lobico_ori, data_feat, all_drugs, drugs, cfold, arch, slr, keepprob, log_file, pred_file):
    print(drugs)
    print(cfold)
    tf.reset_default_graph()
    
    train_drugs = list(set(all_drugs) - set(drugs))
    
    train_index = list(lobico.loc[lobico.Compound.isin(train_drugs)].index)
    valid_index = list(lobico.loc[lobico.Compound.isin(drugs)].index)
    valid_cells = list(lobico.loc[lobico.Compound.isin(drugs)].cell_name)
    valid_drugs = list(lobico.loc[lobico.Compound.isin(drugs)].Compound)
    print(len(valid_index))

    train_feat   = np.asarray(data_feat.loc[train_index])
    valid_feat   = np.asarray(data_feat.loc[valid_index])
    
    train_target = np.asarray(lobico.loc[train_index].value)
    valid_target = np.asarray(lobico.loc[valid_index].value)

    # Standardize initial inputs to 0-1 range
    train_feat, valid_feat = scale_0_1_multiple(train_feat, valid_feat)

    print("SHAPES!!!", train_feat.shape, valid_feat.shape)

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

        main_dict["h0d"] = x

        for l in xrange(len(layers)-1):
            main_dict["W%s"%(l+1)]  = weight_variable([ layers[l] , layers[l+1]])

            main_dict["h%s"%(l+1)]  = dense_lrelu_batch(tf.matmul(main_dict["h%sd"%l], main_dict["W%s"%(l+1)]), 
                                                        phase_train, scope="h%s%s"%(l+1, cfold))
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
        # train_step    = tf.train.AdamOptimizer(slr).minimize(cost)
        train_step    = tf.contrib.opt.NadamOptimizer(slr).minimize(cost)
    
    sess          = tf.InteractiveSession(config=tf.ConfigProto(log_device_placement=True))
    sess.run(tf.global_variables_initializer())

    # Parameters
    batch_size    = 200

    # Validation checks
    max_training   = 2001 # Maximum number of epochs to train for
    count_batch    = 0 # To get tabs of batch number
    count          = 0 # Number of epochs
    check          = 1000 # Batch size during cost checking
    train_checks   = int(math.ceil(train_feat.shape[0] / float(check) )) #n_samples/check
    valid_checks   = int(math.ceil(valid_feat.shape[0] / float(check) ))
    n_batches      = train_feat.shape[0] / batch_size # n_samples/batch_size

    
    print("n_batches: %s"%n_batches)
    print("train_checks: %s"%train_checks)
    while (count < max_training):
        
        batch_all    = random.sample(xrange(train_feat.shape[0]), batch_size)
        feat_batch   = train_feat[batch_all,]
        target_batch = train_target[batch_all]
        
        train_step.run(feed_dict={x: feat_batch, y_: target_batch,
                                  keep_prob: keepprob, phase_train: True})

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
            
            validation     = [sess.run([cost, prediction],
                                        feed_dict =
                                            {x:valid_feat[i*check : (i+1)*check,],
                                             y_:valid_target[i*check : (i+1)*check,],
                                             keep_prob:1.0, phase_train:False})
                              for i in xrange(valid_checks)]
            valid_cost     = np.mean([i[0] for i in validation])
            
            # Calculate ROC per drug in validation set
            p   = [i[1] for i in validation]
            p   = [i for j in p for i in j]
            roc = pd.DataFrame({"cell_name":valid_cells, "Compound":valid_drugs, 
                               "prediction":p})
            roc      = pd.merge(roc, lobico_ori, on=["Compound", "cell_name"])
            # roc_mean = np.mean([roc_auc_score(j.value, j.prediction) for i,j in roc.groupby("Compound")])
            
            # Write to log
            with open(log_file, "a") as f:
                f.write("\n%s\t%s\t%s\t%s\t%s\t%s\t%s"%(cfold, count, batch_size, keepprob,
                                                            train_cost, valid_cost, slr))
            
            # Write predictions
            roc = roc.assign(CFOLD=cfold, Epoch=count)
            roc = roc[["CFOLD", "Epoch","Compound","cell_name","prediction", "value"]]
            with open(pred_file, "a") as f:
                roc.to_csv(f, header=False, sep="\t", index=False)
            
        else:
            count_batch+=1
    sess.close()
    print("Done modeling pack-cell MLP")

def run():
    # Arguments
    out_folder      = "/tigress/zamalloa/CGP_FILES/RESULTS/"
    exp_arch        = sys.argv[1] # The expression features dae (cgp_encode.py) results
    exp_gfilter     = int(sys.argv[2]) # The expression g_filter used (cgp_encode.py)
    drug_arch       = sys.argv[3] # "8_16" # The drug features dae (drug_dae.py) results
    keepprob        = float(sys.argv[4])
    arch            = sys.argv[5]
    slr             = float(sys.argv[6])
    target_features = bool(sys.argv[7]=="True")
    exp_target      = bool(sys.argv[8]=="True")
    string_th       = int(sys.argv[9]) #Irrelevant if both target_features==False & exp_target==False

    # Load data
    lobico, exp_data, drug_data, lobico_ori = load_data(exp_arch, exp_gfilter, drug_arch)

    # Build feature space
    lobico, data_feat   = parse_features(lobico, exp_data, drug_data, target_features, exp_target, string_th)

    # Prep output files (for all folds)
    log_file, pred_file = prep_output_files(out_folder, exp_arch, exp_gfilter, drug_arch, 
                                            keepprob, arch, slr, target_features, exp_target, string_th)

    # Split into folds
    all_drugs, drug_chunks = chunk_split(lobico, folds=5)

    # Train network
    for cfold in xrange(len(drug_chunks)):
        mlp(lobico, lobico_ori, data_feat, all_drugs, drug_chunks[cfold], (cfold+1), 
            arch, slr, keepprob, log_file, pred_file)

    print("DONE")

if __name__ == "__main__":
    run()
