#COUNT OF COMMON GENES COVERED BY EXPRESSION DATA ACROSS PATIENTS
#NEED:
#   Level 3 cancer data
#01/21/14

import subprocess

#Load folder
CANCER_FOLDER=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/Level_3").splitlines()

#Treat files and extract data
PATIENT_DICT={}
for record in CANCER_FOLDER:
    FILE_IN=open("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/Level_3/"+record)
    RECORD=[X.split("\t") for X in FILE_IN.read().splitlines()[2:]]
    FILE_IN.close()
    
    #Filter out invalid expression values
    RECORD=filter(lambda x : x[1]!="null", RECORD)
    
    #Store values in dictionary per patient - patient:gene:exp
    PATIENT_DICT[record]=dict((Y[0], float(Y[1])) for Y in RECORD)

#Pairwise comparisson and correlation, then to file
FILE_OUT1=open("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/PATIENT_CORR_GENE_COUNT", "w")
FILE_OUT1.write("PATIENT_1"+"\t"+"PATIENT_2"+"\t"+"COMMON_GENES")
TOTAL_PATIENTS=len(PATIENT_DICT.keys())

for patient1 in PATIENT_DICT.keys():
    #print (TOTAL_PATIENTS-len(PATIENT_DICT.keys())) *100/float(TOTAL_PATIENTS)
    
    for patient2 in PATIENT_DICT.keys():
        
        #Compare correlation among common genes only
        COMMON_GENES=set(PATIENT_DICT[patient1].keys()) & set(PATIENT_DICT[patient2].keys())
        
        #Store
        FILE_OUT1.write("\n"+patient1+"\t"+patient2+"\t"+str(len(COMMON_GENES)))
    
    #Remove key so we don't recalculate twice
    #del PATIENT_DICT[patient1]

FILE_OUT1.close()