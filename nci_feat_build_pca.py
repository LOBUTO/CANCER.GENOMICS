#nci60_feat_build_pca.py
#Function to build pca training, validation and testing table from nci60 data and pca pre-built model

import cPickle
import random
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale

######################################################################################################
#LOAD INPUTS
n_pcas = int(sys.argv[1])
print(n_pcas)
splits = 0.8

######################################################################################################
#LOAD FILES
PCA = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA/"
PCA_MODEL = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA.MODELS/"
TABLES = "/tigress/zamalloa/TABLES/TCGA.TRAINING/"

with open(PCA_MODEL + "nci60.pca.3.random208.pc1000.pkl", "rb") as md:
    pca_model = cPickle.load(md)

df_test  = pd.read_csv(TABLES + "nci60_all_feat_table_scaled_round.3.csv", sep="\t", nrows=100)
float_cols = [c for c in df_test] # if df_test[c].dtype=="float64"]
float32_cols = {c:np.float32 for c in float_cols}
df  = pd.read_csv(TABLES  + "nci60_all_feat_table_scaled_round.3.csv", sep="\t", engine="c", dtype=float32_cols)
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

#Apply rotation
nci60_pca = np.dot(df, rotation)

#Split as matrix
all_rows = xrange(nci60_pca.shape[0])
train_rows = random.sample(all_rows, int(len(all_rows) * splits))
valid_rows = list(set(all_rows) - set(train_rows))

train_table = nci60_pca[train_rows]
valid_table = nci60_pca[valid_rows]

train_labels = df_labels.iloc[train_rows].as_matrix() #To remove pre-indexes
valid_labels = df_labels.iloc[valid_rows].as_matrix() #To remove pre-indexes

#Apply scaling and labels
std_scale = StandardScaler().fit(train_table)
train_scaled = std_scale.transform(train_table)
valid_scaled = std_scale.transform(valid_table)

train_scaled = pd.DataFrame(train_scaled)
valid_scaled = pd.DataFrame(valid_scaled)

train_labels = pd.DataFrame(train_labels)
valid_labels = pd.DataFrame(valid_labels)
train_labels.columns = ["PCA"]
valid_labels.columns = ["PCA"]

print((train_scaled.shape, valid_scaled.shape))
print((train_labels.shape, valid_labels.shape))

train_scaled = pd.concat([train_labels, train_scaled], axis=1)
valid_scaled = pd.concat([valid_labels, valid_scaled], axis=1)

#Obtain tcga features as PC and scale
tcga_feat = [c for c in df_test]
tcga_feat = tcga[tcga_feat]

tcga_pca = scale(tcga_feat)
tcga_pca = np.dot(tcga_pca, rotation)
tcga_pca = std_scale.transform(tcga_pca) #StandardScaler
tcga_pca = pd.DataFrame(tcga_pca)

tcga_labels = scale(tcga.LIVED.astype(float)) #May get dtype warning
tcga_labels = pd.DataFrame(tcga_labels)
tcga_labels.columns = ["LIVED"]

tcga_pca = pd.concat([tcga_labels, tcga_pca], axis=1)
print(list(tcga_pca.columns.values)[:10])

print("Done executing")
#####################################################################################################
#Pickle files
print(train_scaled.shape)
print(valid_scaled.shape)

with open(TABLES + "nci60_train_matrix" + str(n_pcas) + ".pkl", "wb") as tr:
    cPickle.dump(train_scaled.as_matrix(), tr)

with open(TABLES + "nci60_valid_matrix" + str(n_pcas) + ".pkl", "wb") as vd:
    cPickle.dump(valid_scaled.as_matrix(), vd)

with open(TABLES + "tcga_test_matrix" + str(n_pcas) + ".pkl", "wb") as tc:
    cPickle.dump(tcga_pca.as_matrix(), tc)

with open(TABLES + "nci60_train_scaling" + str(n_pcas) + ".pkl", "wb") as sc:
    cPickle.dump(std_scale, sc)

print("DONE!")
