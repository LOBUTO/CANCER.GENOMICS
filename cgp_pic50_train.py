#cgp_pic50_train.py
import cPickle
import random
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale
import gc

########################################################################################################
#LOAD FILES
CGP_FILES = "/home/zamalloa/Documents/FOLDER/CGP_FILES/"
TRAIN_TABLES = "/home/zamalloa/Documents/FOLDER/TABLES/CGP.TRAINING/"

CGP_FILES = "/tigress/zamalloa/CGP_FILES/"
TRAIN_TABLES = "/tigress/zamalloa/CGP_FILES/CGP_TRAIN_TABLES/"

main_table = pd.read_csv(CGP_FILES + "cgp_cor_pIC50.csv", sep="\t")
nci_table = pd.read_csv(CGP_FILES + "nci60_cgpfeat_cancertable.csv", sep="\t")
nci_stats   = pd.read_csv(CGP_FILES + "nci60_stats.csv", sep="\t")
nci_cgp_dict = pd.read_csv(CGP_FILES + "cgp_to_nci_cell_exp", sep="\t")

target_drug = sys.argv[1] #For testing
nci_boost = bool(sys.argv[2]) #True/False
target_filter = bool(sys.argv[3]) #True/False

#target_cancers = ["Non-Small Cell Lung"] #EXAMPLE for cells related to Erlotinib (nscl)
#target_cancers = ["Leukemia"] #EXAMPLE for cells related to Bosutinib
#target_cancers = ["Leukemia", "Non-Small Cell Lung", "Central Nervous System"] #EXAMPLE for cells related to A-443654
splits = 0.8

print("Done loading files")
########################################################################################################
#EXECUTE
test_table = main_table[main_table.Compound==target_drug]
target_cells = list(set(test_table.cell_name))

main_table = main_table[main_table.Compound!=target_drug]

if target_filter==True:
    main_table = main_table[main_table.cell_name.isin(target_cells)]

if nci_boost==True:

  #nci_cells = list(nci_stats[nci_stats.cancer.isin(target_cancers)]["cell.name"])
  nci_cells = nci_cgp_dict[nci_cgp_dict.Compound==target_drug]["cell_name"]
  nci_table = nci_table[nci_table.cell_name.isin(nci_cells)]

  main_table = main_table.append(nci_table, ignore_index=True)
  main_table.reset_index(drop=True)

gc.collect()
print("Done combining tables")

#Split tables
all_rows = xrange(main_table.shape[0])
train_rows = random.sample(all_rows, int(len(all_rows) * splits))
valid_rows = list(set(all_rows) - set(train_rows))

main_table = main_table.iloc[:,2:].as_matrix()
train_table = main_table[train_rows]
valid_table = main_table[valid_rows]

test_table = test_table.iloc[:,2:].as_matrix()

gc.collect()
print("Done splitting tables")

#Apply scaling
train_labels = train_table[:,0:1]
valid_labels = valid_table[:,0:1]
test_labels  = test_table[:,0:1]

std_scale = StandardScaler().fit(train_table[:,1:])
train_scaled = std_scale.transform(train_table[:,1:])
valid_scaled = std_scale.transform(valid_table[:,1:])
test_scaled = std_scale.transform(test_table[:,1:])

train_scaled = np.concatenate((train_labels, train_scaled), axis=1)
valid_scaled = np.concatenate((valid_labels, valid_scaled), axis=1)
test_scaled  = np.concatenate((test_labels, test_scaled), axis=1)

gc.collect()

print("Done executing")
#####################################################################################################
#Pickle files
print(train_scaled.shape)
print(valid_scaled.shape)
print(test_scaled.shape)

with open(TRAIN_TABLES + "TRAIN." + target_drug + ".pIC50.pkl", "wb") as tr:
    cPickle.dump(train_scaled, tr)

with open(TRAIN_TABLES + "VALID." + target_drug + ".pIC50.pkl", "wb") as vd:
    cPickle.dump(valid_scaled, vd)

with open(TRAIN_TABLES + "TEST." + target_drug + ".pIC50.pkl", "wb") as tc:
    cPickle.dump(test_scaled, tc)

print("DONE!")
