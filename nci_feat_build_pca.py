#nci60_feat_build_pca.py
#Function to build pca training, validation and testing table from nci60 data and pca pre-built model

import cPickle
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

######################################################################################################
#LOAD FILES
PCA = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA/"
PCA_MODEL = "/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA.MODELS/"
TABLES = "/tigress/zamalloa/TABLES/TCGA.TRAINING/"

with open(PCA_MODEL + "nci60.pca.random134.pc200.pkl") as md:
    pca_model = cPickle.load(md)

df_test  = pd.read_csv(TABLES + "nci60_all_feat_table_scaled_round.csv", sep="\t", nrows=100)
float_cols = [c for c in df_test if df_test[c].dtype=="float64"]
float32_cols = {c:np.float32 for c in float_cols}
df  = pd.read_csv(TABLES  + "nci60_all_feat_table_scaled_round.csv", sep="\t", engine="c", dtype=float32_cols)

tcga = pd.read_csv(TABLES + "062816_tcga_all_feat_table.csv", sep="\t")

#####################################################################################################
#EXECUTE

#Obtain nci60 features as PC
nci60_pca = pca_model.transform(df)

#Obtain tcga features as PC
tcga_feat = [c for c in df]
tcga_feat = tcga[tcga_feat]

tcga_pca = scale(tcga_feat)
tcga_pca = pca_model.transform(tcga_pca)

tcga_labels = tcga.LIVED

#Write to table
