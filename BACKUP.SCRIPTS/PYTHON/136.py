#EXTRACT DISTANCES BASED ON DIRECT CORRELATION
#RE-DO OF 134.py 1-corr used as distance instead of hierarchical distance
#020414
import numpy as np
import pickle

FOLDER1="DATABASES/CANCER_DATA/TCGA/"

#Load TOP BRCA PATIENTS - CAN MODIFY WHICH ARE TOP PATIENTS BASED ON HEIGHT
FILE_IN1=open(FOLDER1+"FROM_R_ANALYSIS/012514_BRCA_TOP_CORR_PATIENTS_490")
TOP_PATIENTS=[x for x in FILE_IN1.read().splitlines()]
FILE_IN1.close()

#Load patient pairwise correlation
FILE_IN2=open(FOLDER1+"EXP_GENE/BRCA/022414_PATIENT_CORR")
PATIENT_COR=[y.split("\t") for y in FILE_IN2.read().splitlines()]
FILE_IN2.close()

#Create $ correlation dictionary
PATIENT_COR_DICT=dict((Z[0]+"$"+Z[1], float(Z[2])) for Z in PATIENT_COR)
PATIENT_COR_DICT.update( dict((Z[1]+"$"+Z[0], float(Z[2])) for Z in PATIENT_COR) )

ALL_PATIENTS=list(set([X[0] for X in PATIENT_COR]))

#Load ID to SAMPLE data of the form - dict patient_id:sample_id - STRING HAS EXTRA CHARACTERS!!!! - THIRD
FILE_IN3=open("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/9078bd08-7662-45b9-ab3b-25ddf04f2897/FILE_SAMPLE_MAP.txt")
PATIENT_TO_SAMPLE_DICT=dict((Z.split("\t")[0], Z.split("\t")[1][:16]) for Z in FILE_IN3.read().splitlines()[1:])
FILE_IN3.close()

#GET TO PICKLE FORMAT OF OUT_FILE IN 134 and repeat 135 with it

#Create pairwise distances based on correlation to core
PICKLE_OUT1=open("DATABASES/CANCER_DATA/TCGA/ALL_PATIENT_COR_DIST_DICT_1_0.pi", "w")
DICT_OUT1={}
for patient in ALL_PATIENTS:
    print patient
    
    if patient in TOP_PATIENTS:
        DICT_OUT1[PATIENT_TO_SAMPLE_DICT[patient]]=0.0
    
    else:
        DICT_OUT1[PATIENT_TO_SAMPLE_DICT[patient]]=1.0-np.mean([PATIENT_COR_DICT[patient+"$"+top] for top in TOP_PATIENTS])

pickle.dump(DICT_OUT1,PICKLE_OUT1)
PICKLE_OUT1.close()
             
