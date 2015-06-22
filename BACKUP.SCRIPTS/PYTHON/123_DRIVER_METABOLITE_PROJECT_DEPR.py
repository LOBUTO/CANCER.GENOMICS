#DRIVER METABOLITE PROJECT - Script not in use, modified what sifnificant means in 12/19/13
#12/18/13
#From R - Significant metabolites per cancer
#Process to find how many patients per cancer have at least one of them

import pickle

#Load and process Significant metabolites file
FILE_IN1=open("NETWORK/FROM_R_ANALYSIS/121813_SIG_MET_PER_CANCER_100")
SIG_PER_CANCER=[X.split("\t")[:2] for X in FILE_IN1.read().splitlines()]
FILE_IN1.close()

SIG_CANCER_DICT=dict((record[1],[]) for record in SIG_PER_CANCER)

for record in SIG_PER_CANCER: 
    SIG_CANCER_DICT[record[1]]=SIG_CANCER_DICT[record[1]] + [record[0]] #Of the form {cancer:[SIGNIFICANT_METABOLITES]} 

#Load PMV dict of the form {cancer:patient:metabolite:normalized count}
PICKLE_IN1=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PMV_NETWORK100.pi")
DICT_PMV=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

for cancer in DICT_PMV.keys():
    
    for patient in DICT_PMV[cancer].keys():
        
        if len(set(DICT_PMV[cancer][patient].keys()) & set(SIG_CANCER_DICT[cancer]))==0:
            del DICT_PMV[cancer][patient]
    
    print cancer, len(DICT_PMV[cancer].keys())
            
            
