#nci60_feat_build_pca.py
#Function to build pca training, validation and testing table from nci60 data and pca pre-built model

import cPickle
import random
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

######################################################################################################
#LOAD INPUTS
n_pcas = 300
splits = 0.8

######################################################################################################
#LOAD FILES
PCA = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA/"
PCA_MODEL = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA.MODELS/"
TABLES = "/tigress/zamalloa/TABLES/TCGA.TRAINING/"

with open(PCA_MODEL + "nci60.pca.random214.pc500.pkl", "rb") as md:
    pca_model = cPickle.load(md)

df_test  = pd.read_csv(TABLES + "nci60_all_feat_table_scaled_round.2.csv", sep="\t", nrows=100)
float_cols = [c for c in df_test] # if df_test[c].dtype=="float64"]
float32_cols = {c:np.float32 for c in float_cols}
df  = pd.read_csv(TABLES  + "nci60_all_feat_table_scaled_round.2.csv", sep="\t", engine="c", dtype=float32_cols)
print(df.shape)

df_labels = pd.read_csv(PCA + "nci60.pca.labels.all.csv", sep="\t")
df_labels.columns = ["PCA"]
print(df_labels.shape)

tcga = pd.read_csv(TABLES + "062816_tcga_all_feat_table.csv", sep="\t")
print(tcga.shape)

print("Done loading files")
#####################################################################################################
#EXECUTE

#Obtain rotation matrix (Eigenvectors)
rotation = pca_model.components_[:n_pcas]
rotation = rotation.transpose()

#Obtain nci60 features as PC and scale
nci60_pca = np.dot(df, rotation)
nci60_pca = scale(nci60_pca)
nci60_pca = pd.DataFrame(nci60_pca)

nci60_pca = pd.concat([df_labels, nci60_pca], axis=1)
print(list(nci60_pca.columns.values)[:10])
with open(TABLES + "nci_pca.pkl", "wb") as pca:
    cPickle.dump(nci60_pca, pca) #Store for now

#Obtain tcga features as PC and scale
tcga_feat = [c for c in df_test]
tcga_feat = tcga[tcga_feat]

tcga_pca = scale(tcga_feat)
tcga_pca = np.dot(tcga_pca, rotation)
tcga_pca = scale(tcga_pca)
tcga_pca = pd.DataFrame(tcga_pca)

tcga_labels = scale(tcga.LIVED.astype(float)) #May get dtype warning
tcga_labels = pd.DataFrame(tcga_labels)
tcga_labels.columns = ["LIVED"]

tcga_pca = pd.concat([tcga_labels, tcga_pca], axis=1)
print(list(tcga_pca.columns.values)[:10])

print("Done executing")
#####################################################################################################
#Pickle files
all_rows = xrange(nci60_pca.shape[0])
train_rows = random.sample(all_rows, int(len(all_rows) * splits))
valid_rows = list(set(all_rows) - set(train_rows))

train_table = nci60_pca.iloc[train_rows]
valid_table = nci60_pca.iloc[valid_rows]
print(train_table.shape)
print(valid_table.shape)

with open(TABLES + "nci60_train_matrix.pkl", "wb") as tr:
    cPickle.dump(train_table.as_matrix(), tr)

with open(TABLES + "nci60_valid_matrx.pkl", "wb") as vd:
    cPickle.dump(valid_table.as_matrix(), vd)

with open(TABLES + "tcga_test_matrix.pkl", "wb") as tc:
    cPickle.dump(tcga_pca.as_matrix(), tc)

print("DONE!")
