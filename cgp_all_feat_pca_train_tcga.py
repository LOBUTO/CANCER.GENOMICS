#cgp_all_feat_pca_train_tcga.py
#Function to prep train, valid and target pca tables for training model, this will target for tcga data

#####################################################################################################################
import cPickle
import random
import sys
import numpy as np
import pandas as pd
import gc
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

######################################################################################################
#LOAD INPUTS
n_pcas = int(sys.argv[1])
print(n_pcas)
splits = 0.8

######################################################################################################
#LOAD FILES
FOLDER = "/tigress/zamalloa/CGP_FILES/"
TRAIN_TABLES = "/tigress/zamalloa/CGP_FILES/CGP_TRAIN_TABLES/"
MODELS = "/tigress/zamalloa/CGP_FILES/CGP_MODELS/"

with open(FOLDER + "cgp_all_feat", "r") as cg:
    cgp_table = pd.read_csv(cg, sep = "\t")

tcga = pd.read_csv("/tigress/zamalloa/TABLES/TCGA.TRAINING/" + "062816_tcga_all_feat_table.csv", sep="\t")

######################################################################################################
#EXECUTE

#Split tables
all_feat = list(cgp_table.columns.values[3:])
msk = np.random.rand(len(cgp_table)) < splits

train_table = cgp_table[msk]
valid_table = cgp_table[~msk]

tcga_feat = tcga[all_feat]
tcga_feat = scale(tcga_feat) #PRE-scaling, just in case

del cgp_table #To save memory
gc.collect()
print(train_table.shape)
print(valid_table.shape)

#Calculate PCA on training data and store model (If not previously provided)
if len(sys.argv)>2 :

    with open(sys.argv[2], "rb") as ff:
        pca = cPickle.load(ff)

else:

    pca = PCA(n_components=1000)
    pca.fit(train_table.iloc[:,3:])

    with open(MODELS + "cgp_pca_1000.pkl", "wb") as ff:
        cPickle.dump(pca, ff)

var1=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
print(var1)

#Apply rotation
rotation = pca.components_[:n_pcas]
rotation = rotation.transpose()

train_labels = pd.DataFrame({"NORM_AUC" : list(train_table["NORM.AUC"])})

train_table = scale(np.dot(train_table.iloc[:,3:], rotation)) #Post scaling
train_table = pd.DataFrame(train_table)
print(train_labels.shape, train_table.shape)
train_labels.reset_index(drop=True, inplace=True); train_table.reset_index(drop=True, inplace=True)
train_table = pd.concat([train_labels, train_table], axis=1)

valid_labels = pd.DataFrame({"NORM_AUC" : list(valid_table["NORM.AUC"])})

valid_table = scale(np.dot(valid_table.iloc[:,3:], rotation)) #Post scaling
valid_table = pd.DataFrame(valid_table)
print(valid_labels.shape, valid_table.shape)
valid_labels.reset_index(drop=True, inplace=True); valid_table.reset_index(drop=True, inplace=True)
valid_table = pd.concat([valid_labels, valid_table], axis=1)

test_labels = pd.DataFrame({"LIVED" : list(scale(tcga.LIVED.astype(float))) })

test_table = scale(np.dot(tcga_feat, rotation)) #Post scaling
test_table = pd.DataFrame(test_table)

print(test_labels.shape, test_table.shape)
print(test_labels.iloc[:5,:])
print(test_table.iloc[:5,:5])
test_labels.reset_index(drop=True, inplace=True) ; test_table.reset_index(drop=True, inplace=True)
test_table = pd.concat([test_labels, test_table], axis=1)

print(train_table.shape)
print(valid_table.shape)
print(test_table.shape)

drug = "TCGA"
with open(TRAIN_TABLES + "cgp_train_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as tr:
    cPickle.dump(train_table.as_matrix(), tr)

with open(TRAIN_TABLES + "cgp_valid_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as vd:
    cPickle.dump(valid_table.as_matrix(), vd)

with open(TRAIN_TABLES + "cgp_test_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as tc:
    cPickle.dump(test_table.as_matrix(), tc)

print("DONE!")
