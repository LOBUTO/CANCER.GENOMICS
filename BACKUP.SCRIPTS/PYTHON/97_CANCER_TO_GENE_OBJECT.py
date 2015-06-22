#Separate cancer genes by cancer type

import subprocess, pickle
from FUNCTIONS import unique

#Get genes in stats
GENES_STATS=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/geneStats").splitlines()

#CANCER TYPES
CANCER_DICT={"LUSC":[], "READ":[], "GBM":[], "KIRC":[], "UCEC":[], "OV":[], "BRCA":[], "COAD":[]}

#Assign to CANCER TYPES
for gene in GENES_STATS:
    FILE_IN=open("DATABASES/CANCER_DATA/geneStats/%s"%gene)
    RECORDS=FILE_IN.read().splitlines()
    FILE_IN.close()
    
    for record in RECORDS:
        CANCER_DICT[record.split()[4].strip()]=CANCER_DICT[record.split()[4].strip()]+[gene.split(".")[0]]

for cancer in CANCER_DICT.iterkeys():
    CANCER_DICT[cancer]=unique(CANCER_DICT[cancer])
    print cancer, CANCER_DICT[cancer]

#Store as object
PICKLE_OUT=open("DATABASES/OBJECTS/DICT_CANCER_TO_GENES.pi", "w")
pickle.dump(CANCER_DICT, PICKLE_OUT)

for cancer in CANCER_DICT.iterkeys():
    print cancer, len(CANCER_DICT[cancer])
    