#BACKGROND MODEL PROJECT 1
#Build matrix of top BRCA patients for PCA analysis in R
#012514 

FOLDER1="DATABASES/CANCER_DATA/TCGA/"

#Load TOP BRCA PATIENTS
FILE_IN1=open(FOLDER1+"FROM_R_ANALYSIS/012514_BRCA_TOP_CORR_PATIENTS")
TOP_PATIENTS=[x for x in FILE_IN1.read().splitlines()]
FILE_IN1.close()

#Load expression data into dict
TOP_PATIENT_DICT={}
for patient in TOP_PATIENTS:
    #Load expression file
    FILE_IN1=open(FOLDER1+"EXP_GENE/BRCA/Level_3/%s"%patient)
    PATIENT_EXP=[X.split("\t") for X in FILE_IN1.read().splitlines()[2:]]
    FILE_IN1.close()
    
    #Add to dictionary
    TOP_PATIENT_DICT[patient]=dict((Y[0], Y[1]) for Y in PATIENT_EXP if Y[1]!="null")
    
#Find common genes across all TOP PATIENTS
COMMON_GENES=list(set.intersection(*[set(TOP_PATIENT_DICT[patients].keys()) for patients in TOP_PATIENT_DICT.keys()]))
print len(COMMON_GENES) #Common genes across TOP PATIENTS are 17670

#Create matrix for PCA and write out to file
FILE_OUT1=open(FOLDER1+"FOR_R_ANALYSIS/012514_BRCA_TOP_PATIENTS_FOR_PCA", "w")
FILE_OUT1.write("PATIENT")
for gene in COMMON_GENES:
    FILE_OUT1.write("\t"+gene)

for patient in TOP_PATIENT_DICT.keys():
    FILE_OUT1.write("\n"+patient)
    
    for gene in COMMON_GENES:
        FILE_OUT1.write("\t"+str(TOP_PATIENT_DICT[patient][gene]))

FILE_OUT1.close()
