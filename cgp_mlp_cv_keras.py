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
from keras.optimizers import Adam,Nadam,SGD,Adagrad,Adadelta,RMSprop
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
from keras.utils import to_categorical
from keras import backend as K
from keras import metrics

# Import functions
sys.path.insert(0, '/home/ubuntu/GIT')

from paper_1 import *

def keras_auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc

def exp_encode_parse(gsea, g_filter, arch, dropkeep):
    in_folder = "/home/ubuntu/GSEA_FILES/RESULTS/"
    
    main_list = []
    for i in ["train", "valid", "test"]:
        in_file = "{}CGP_{}_{}_deepautoencoder_{}_{}_{}.txt".format(in_folder,gsea,g_filter,arch,dropkeep,i)
        
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

def drug_encode_parse(arch, dropkeep):
    in_folder = "/home/ubuntu/CGP_FILES/RESULTS/"
    
    main_list = []
    for i in ["train", "valid", "test"]:

        in_file = "{}{}_deepautoencoder_{}_{}_{}.txt".format(in_folder, "cgp_drugs", arch, dropkeep,i)
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

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
        
        class_weight = {0:1.0, 1:zero_weight}
#         class_weight = {0:zero_weight, 1:1.0}
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

    model.add(Dense(layers[1], activation='relu', input_dim=layers[0]))
    model.add(BatchNormalization())
    model.add(Dropout(keepprob))    

    for l in layers[2:]:
        model.add(Dense(l, activation='relu'))
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

def drug_feat_target_cv_split(target_table, feat_table, drug, cells):

    # Index based on training sample
    all_index    = list(target_table.index)
    valid_index  = list(target_table.loc[(target_table.Compound==drug) & (target_table.cell_name.isin(cells))].index)
    print("valid_index: ", valid_index)
    train_index  = list(set(all_index) - set(valid_index))

    train_feat   = np.asarray(feat_table.loc[train_index])
    valid_feat   = np.asarray(feat_table.loc[valid_index])
    
    train_target = np.asarray(target_table.loc[train_index].value)
    valid_target = np.asarray(target_table.loc[valid_index].value)

    # Return true_values
    true_values  = target_table.loc[valid_index]["value"].values

    return train_feat, valid_feat, train_target, valid_target, true_values

def prep_output_files_keras_cv(out_folder, exp_arch, exp_gfilter, drug_arch, keepprob, arch, slr, 
                                 target_features, exp_target, encoded_target, encoded_target_arch, string_th, gdsc_variants, problem, drug):
    # Prep output files
    pred_file     = "{}cgpmlp_exp_drug_reg_cv_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_txt".format(out_folder, "drugscv", exp_arch, exp_gfilter, 
                                               drug_arch, keepprob, arch, slr, target_features, exp_target, 
                                               encoded_target, encoded_target_arch, string_th,gdsc_variants, problem, drug)
    if problem=="classification":
        with open(pred_file, "w") as f:
            f.write("Compound\tcell_name\ttrue_value\tprob_0\tprob_1\tbest_metric")
    elif problem=="regression":
        with open(pred_file, "w") as f:
            f.write("Compound\tcell_name\ttrue_value\tprediction\tbest_metric")

    return pred_file

def fix_cell_drug_feat_indices(feat_main, exp_feat, drug_feat, target_feat):

    common_cells = set(list(feat_main.cell_name)).intersection(exp_feat.index)
    common_drugs = set(list(feat_main.Compound)).intersection(drug_feat.index)
    common_drugs = set(list(feat_main.Compound)).intersection(target_feat.index)
    feat_main    = feat_main.loc[feat_main.Compound.isin(common_drugs)]

    data_feat    = pd.concat([exp_feat.loc[feat_main.cell_name,].reset_index(drop=True),
                              drug_feat.loc[feat_main.Compound,].reset_index(drop=True),
                              target_feat.loc[feat_main.Compound,].reset_index(drop=True)], axis=1)
    return feat_main, data_feat

def parse_features(lobico, exp_data, mut_data, drug_data, target_features, exp_target, encoded_target, encoded_target_arch, string_th, gdsc_variants):
    # Build feature space
    
    # Build cell line feature space first
    if gdsc_variants==True:
        print("Using gdsc variant features")
        cell_data = pd.concat([exp_data, mut_data.loc[exp_data.index,]], 
                              axis=1)
    else:
        print("Not using gdsc variant features")
        cell_data = exp_data
    print("cell features: ", cell_data.shape[1])
        
    # Do we need to add target features??
    if target_features==True:
        print("Using drug target features")
        pp_hugo             = string_pp_to_hugo("/home/ubuntu/STRING/9606.protein.links.v10.5.txt",
                                                "STRING/9606.protein.aliases.v10.5.txt")
        drug_target         = process_drug_target("/home/ubuntu/CGP_FILES/Screened_Compounds.csv")

        if exp_target==True:
            print("Using target expression features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="table")
            lobico, target_feat = string_target_expression_features(drug_target, "CGP_FILES/070818_cgp_exp.txt", lobico)

            data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                        drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                        target_feat.reset_index(drop=True)], axis=1)
        else:
            print("Using target binary features")
            
            if encoded_target==True:
                print("Using encoded features")
                target_feat       = pd.read_csv("/home/ubuntu/CGP_FILES/drug_target_encoded_binary_{}_{}.txt".format(string_th, encoded_target_arch), 
                                                sep="\t", index_col=0)
                lobico, data_feat = fix_cell_drug_feat_indices(lobico, cell_data, drug_data, target_feat)

            else:
                print("Using non-encoded features")
                drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="pivot")
                lobico, target_feat = target_lobico_parse(lobico, drug_target)

                data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                            drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                            target_feat.loc[lobico.Compound].reset_index(drop=True)], axis=1)
        
    else:
        print("No target features used")
        data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                    drug_data.loc[lobico.Compound,].reset_index(drop=True)], axis=1)
    
    lobico        = lobico.reset_index(drop=True)
    print(data_feat.head())
    print("All data:", data_feat.shape)

    return lobico, data_feat

def load_data(exp_arch, exp_gfilter, drug_arch, problem="regression"):
    # problem: regression, classification
    # Data for ROC
    lobico_ori = process_paper_lobico("/home/ubuntu/PAPERS/LORIO/lobico.csv")
    lobico_ori["Compound"] = [Function_drug_name_zhang_to_gdsc(i, "ztg") for i in lobico_ori["Compound"]]

    ####### Prep data #####
    # Expression features
    exp_data      = exp_encode_parse("c2setcover", exp_gfilter, exp_arch, 0.5)
    
    # Mutation features
    mut_data      = process_gdsc_variant("/home/ubuntu/CGP_FILES/gdsc_WES_variants.csv", output="matrix", th=5)
    
    # Drug features
    if drug_arch=="original":
        drug_data     = pd.read_csv("/home/ubuntu/PAPERS/ZHANG_2018/drug_feats_proc.txt", sep="\t", index_col=0)
        
    elif drug_arch=="pubchem_smiles":
        drug_data     = pd.read_csv("/home/ubuntu/CGP_FILES/pubchem_smiles.txt", index_col=0)
        
    elif drug_arch=="original_pubchem_smiles":
        print("using original and pubchem smiles")
        original  = pd.read_csv("/home/ubuntu/PAPERS/ZHANG_2018/drug_feats_proc.txt", sep="\t", index_col=0)
        pubchem   = pd.read_csv("/home/ubuntu/CGP_FILES/pubchem_smiles.txt", index_col=0)
        drug_data = pd.merge(original, pubchem, left_index=True, right_index=True)
        
    elif "_arch_" in drug_arch:
        print("using combined datasets")
        c_list    = [i for i in drug_arch.split("_arch_")]
        pubchem   = pd.read_csv("/home/ubuntu/CGP_FILES/pubchem_smiles.txt", index_col=0)
        encoded   = drug_encode_parse(c_list[1], 0.5)
        drug_data = pd.merge(pubchem, encoded, left_index=True, right_index=True)
    else:
        drug_data     = drug_encode_parse(drug_arch, 0.5)

    # Is this for classification or regression? 
    if problem=="regression":
        lobico = cgp_act_post_process(pd.read_csv("/home/ubuntu/CGP_FILES/v17.3_fitted_dose_response.csv"),
                                        zscoring=False)
    elif problem=="classification":
        lobico = lobico_ori

    lobico, exp_data, mut_data, drug_data = drug_exp_lobico_parse_v2(lobico, exp_data, mut_data, drug_data, "ztg")

    return lobico, exp_data, mut_data, drug_data, lobico_ori

def run(exp_arch, exp_gfilter, drug_arch, keepprob, arch, slr, target_features, exp_target,
        encoded_target, encoded_target_arch,
        string_th, gdsc_variants, problem, batch_size, epochs, drug):
    
    out_folder          = "/home/ubuntu/CGP_FILES/RESULTS/"
     
    # Load data
    lobico, exp_data, mut_data, drug_data, lobico_ori = load_data(exp_arch, exp_gfilter, drug_arch, problem)

    # Build feature space
    lobico, data_feat   = parse_features(lobico, exp_data, mut_data, drug_data, target_features, exp_target, 
                                        encoded_target, encoded_target_arch, string_th, gdsc_variants)

    # Prep output files (for all folds)
    pred_file             = prep_output_files_keras_cv(out_folder, exp_arch, exp_gfilter, drug_arch, 
                                            keepprob, arch, slr, 
                                            target_features, exp_target, encoded_target, encoded_target_arch, string_th, 
                                            gdsc_variants, problem, drug)

    all_cells           = list(set(lobico.loc[lobico.Compound==drug].cell_name))
    cell_chunks         = [list(i) for i in np.array_split(all_cells, 10)]
    
    count = 0
    for cells in cell_chunks:
        count+=1
        print("Sample: ", cells, count/10.0)
        # Split into training and validation
        x_train , x_valid , y_train , y_valid, true_values = drug_feat_target_cv_split(lobico, data_feat, drug, cells)

        # Run model per cell
        results, model, early_stopping = keras_mlp(x_train , x_valid , y_train , y_valid,
                                                   arch, slr, keepprob, problem,
                                                   pred_file, batch_size, epochs)
        best_metric = np.min(results.history["val_loss"]) #Given that accuracy is the early stopping metric

        # Predict
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
    drug                = "Sunitinib"
    exp_arch            = "2" # The expression features dae (cgp_encode.py) results
    exp_gfilter         = int(10) # The expression g_filter used (cgp_encode.py)
    drug_arch           = "8_16" # "pubchem_arch_8_16" #"8_16" 
    arch                = "1_4_4"
    keepprob            = float(0.5)
    slr                 = 0.0001 #0.0000001
    target_features     = False
    exp_target          = False #If target features true, then we will use regression features, otherwise binary features
    encoded_target      = False
    encoded_target_arch = "2_16" #Only used if encoded_target==True
    string_th           = 550 #Irrelevant if both target_features==False & exp_target==False
    gdsc_variants       = False
    problem             = "classification" #regression
    batch_size          = 200
    epochs              = 20
    
    run(exp_arch, exp_gfilter, drug_arch, keepprob,
        arch, slr, target_features, exp_target, encoded_target, encoded_target_arch,
        string_th, gdsc_variants, problem, batch_size, epochs,
        drug)
    print("DONE")

if __name__ == "__main__":
    run_wrapper()