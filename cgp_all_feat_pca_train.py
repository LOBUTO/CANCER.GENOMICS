#cgp_all_feat_pca_train.py
#Function to prep train, valid and target pca tables for training model per drug

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
drug="Erlotinib"

######################################################################################################
#LOAD FILES
FOLDER = "/tigress/zamalloa/CGP_FILES/"
TRAIN_TABLES = "/tigress/zamalloa/CGP_FILES/CGP_TRAIN_TABLES/"
MODELS = "/tigress/zamalloa/CGP_FILES/CGP_MODELS/"

# df_test  = pd.read_csv(FOLDER + "cgp_all_feat", sep="\t", nrows=100)
# float_cols = [c for c in df_test if df_test[c].dtype=="float64"]
# float32_cols = {c:np.float32 for c in float_cols}
# df  = pd.read_csv(FOLDER + "cgp_all_feat", sep="\t", engine="c", dtype=float32_cols)

with open(FOLDER + "cgp_all_feat", "r") as cg:
    cgp_table = pd.read_csv(cg, sep = "\t")

######################################################################################################
#EXECUTE

#Split tables
test_table = cgp_table[cgp_table["DRUG"]==drug]
cgp_table = cgp_table[cgp_table["DRUG"]!=drug]
msk = np.random.rand(len(cgp_table)) < splits

train_table = cgp_table[msk]
valid_table = cgp_table[~msk]

del cgp_table #To save memory
gc.collect()
print(train_table.shape)
print(valid_table.shape)

#Calculate PCA on training data and store model
pca = PCA(n_components=1000)
pca.fit(train_table.iloc[:,3:])
var1=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
print(var1)

with open(MODELS + "cgp_pca_1000.pkl", "wb") as ff:
    cPickle.dump(pca, ff)

#Apply rotation
rotation = pca.components_[:n_pcas]
rotation = rotation.transpose()

train_table = np.dot(train_table[:,3:], rotation)
train_table = scale(train_table)
train_table = pd.DataFrame(train_table)
train_table = pd.concat([train_table.iloc[:,:3], train_table], axis=1)

valid_table = np.dot(valid_table[:,3:], rotation)
valid_table = scale(valid_table)
valid_table = pd.DataFrame(valid_table)
valid_table = pd.concat([valid_table.iloc[:,:3], valid_table], axis=1)

test_table = np.dot(test_table[:,3:], rotation)
test_table = scale(test_table)
test_table = pd.DataFrame(test_table)
test_table = pd.concat([test_table.iloc[:,:3], test_table], axis=1)

print(train_table.shape)
print(valid_table.shape)
print(test_table.shape)

with open(TRAIN_TABLES + "cgp_train_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as tr:
    cPickle.dump(train_table.as_matrix(), tr)

with open(TRAIN_TABLES + "cgp_valid_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as vd:
    cPickle.dump(valid_table.as_matrix(), vd)

with open(TRAIN_TABLES + "cgp_test_matrix_" + drug + str(n_pcas) + ".pkl", "wb") as tc:
    cPickle.dump(test_table.as_matrix(), tc)

print("DONE!")
