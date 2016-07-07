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

train_labels = train_table["NORM.AUC"]
train_labels = pd.DataFrame(train_labels)
train_labels.columns = ["NORM_AUC"]
train_table = np.dot(train_table.iloc[:,3:], rotation)
train_table = scale(train_table)
train_table = pd.DataFrame(train_table)
print(train_labels.shape, train_table.shape)
train_labels.reset_index(drop=True, inplace=True)
train_table.reset_index(drop=True, inplace=True)
train_table = pd.concat([train_labels, train_table], axis=1)

valid_labels = valid_table["NORM.AUC"]
valid_labels = pd.DataFrame(valid_labels)
valid_labels.columns = ["NORM_AUC"]
valid_table = np.dot(valid_table.iloc[:,3:], rotation)
valid_table = scale(valid_table)
valid_table = pd.DataFrame(valid_table)
print(valid_labels.shape, valid_table.shape)
valid_labels.reset_index(drop=True, inplace=True)
valid_table.reset_index(drop=True, inplace=True)
valid_table = pd.concat([valid_labels, valid_table], axis=1)

test_labels = test_table["NORM.AUC"]
test_labels = pd.DataFrame(test_labels)
test_labels.columns = ["NORM_AUC"]
test_table = np.dot(test_table.iloc[:,3:], rotation)
test_table = scale(test_table)
test_table = pd.DataFrame(test_table)
print(test_labels.shape, test_table.shape)
print(test_labels.iloc[:5,:])
print(test_table.iloc[:5,:5])
test_labels.reset_index(drop=True, inplace=True)
test_table.reset_index(drop=True, inplace=True)
test_table = pd.concat([test_labels, test_table], axis=1)

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
