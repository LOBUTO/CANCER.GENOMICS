#BACKGROND MODEL PROJECT 3
#TOP and NON-TOP patient_ids to sample_ids with distance to TOP
#012614

import pickle

#Load patient numeric dictionary from R - dict of the form #:patient_id - SECOND
FILE_IN1=open("DATABASES/CANCER_DATA/TCGA/FROM_R_ANALYSIS/012614_BRCA_PATIENT_DICT_NUMERIC")
PATIENT_NUM_DICT=dict((X.split("\t")[0], X.split("\t")[1]) for X in FILE_IN1.read().splitlines())
FILE_IN1.close()

#Load NON-TOP distances of the form - dict numeric_node:distance_to_top - FIRST
FILE_IN2=open("DATABASES/CANCER_DATA/TCGA/012614_BRCA_NON_TOP_NUMERIC_DISTANCE_TO_TOP")
NON_TOP_DIST_DICT=dict((Y.split("\t")[0], Y.split("\t")[1]) for Y in FILE_IN2.read().splitlines())
FILE_IN2.close()

#Load ID to SAMPLE data of the form - dict patient_id:sample_id - STRING HAS EXTRA CHARACTERS!!!! - THIRD
FILE_IN3=open("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/9078bd08-7662-45b9-ab3b-25ddf04f2897/FILE_SAMPLE_MAP.txt")
PATIENT_TO_SAMPLE_DICT=dict((Z.split("\t")[0], Z.split("\t")[1][:16]) for Z in FILE_IN3.read().splitlines()[1:])
FILE_IN3.close()

#Load TOP_PATIENTS patient_id - FOURTH
FILE_IN4=open("DATABASES/CANCER_DATA/TCGA/FROM_R_ANALYSIS/012514_BRCA_TOP_CORR_PATIENTS")
TOP_PATIENTS=[W for W in FILE_IN4.read().splitlines()]
FILE_IN4.close()

#Write out file of the form - dict sample_id:distance
DICT_OUT1={}
PICKLE_OUT1=open("DATABASES/CANCER_DATA/TCGA/ALL_PATIENT_DIST_DICT.pi", "w")

#First non_top
for patient in NON_TOP_DIST_DICT.keys():
    
    PATIENT_ID=PATIENT_NUM_DICT[patient]
    SAMPLE_ID=PATIENT_TO_SAMPLE_DICT[PATIENT_ID]
    DICT_OUT1[SAMPLE_ID]=NON_TOP_DIST_DICT[patient]

#Then top
for patients in TOP_PATIENTS:
    
    SAMPLE_IDS=PATIENT_TO_SAMPLE_DICT[patients]
    DICT_OUT1[SAMPLE_IDS]=0

pickle.dump(DICT_OUT1, PICKLE_OUT1)
PICKLE_OUT1.close()

