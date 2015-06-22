#DISTANCE AND MUTATION INFO OF EXPRESSION PATIENTS
#RE-DO of 135.py using info obtained from 136.py
#020414

import pickle
 
#Load dict - patient:distance - BASED ON 1-CORR
PICKLE_IN1=open("DATABASES/CANCER_DATA/TCGA/ALL_PATIENT_COR_DIST_DICT_1_0.pi")
PATIENT_DIST_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

#Load mutation data
FILE_IN1=open("DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/012714_BRCA_sample.mut")
BRCA_MUT=[X.split("\t") for X in FILE_IN1.read().splitlines()]
FILE_IN1.close()
"""
FILE_IN2=open("DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/BRCA.mut")
BRCA_MUT=BRCA_MUT+[X.split("\t") for X in FILE_IN2.read().splitlines()]
FILE_IN2.close()
"""
BRCA_MUT=list(set(["_".join(X) for X in BRCA_MUT]))
BRCA_MUT=[X.split("_") for X in BRCA_MUT]


#Do patient:#genes, patient:#isoforms dict, patient:#total mutations (not double counting isoforms), and patient:#mutations of first isoform
PATIENTS=list(set([Y[1].strip() for Y in BRCA_MUT]))
PATIENT_GENE_DICT={}
PATIENT_ISO_DICT={}
PATIENT_MUT_DICT={}

PATIENT_FIRSTISO_DICT={}
BRCA_MUT_FIRSTISO=filter(lambda x: ".001" in x[0] ,BRCA_MUT)

print len(set(PATIENTS) & set(PATIENT_DIST_DICT.keys())) 

for patient in PATIENTS:
    PATIENT_GENE_DICT[patient]= len(set([Y[0].split(".")[0] for Y in BRCA_MUT if Y[1]==patient]))
    PATIENT_ISO_DICT[patient]= len(set([Z[0] for Z in BRCA_MUT if Z[1]==patient]))
    
    PATIENT_MUT_DICT[patient]= len(set([ "_".join( [W[0].split(".")[0], W[2],W[3],W[4]] ) for W in BRCA_MUT if W[1]==patient ]))
    
    PATIENT_FIRSTISO_DICT[patient]= len(filter(lambda Q: Q[1]==patient, BRCA_MUT_FIRSTISO))

#Write out to table for R analysis - CAN MODIFY BASED ON ORIGINAL HIERARCHICAL HEIGHT ON 136.py
FILE_OUT1=open("DATABASES/CANCER_DATA/TCGA/FOR_R_ANALYSIS/012614_BRCA_PATIENT_COR_DIST_MUTATIONS_1_0","w")
FILE_OUT1.write("PATIENT"+"\t"+"DISTANCE"+"\t"+"GENES"+"\t"+"ISOFORMS"+"\t"+"ALL_MUTATIONS"+"\t"+"FIRST_ISOFORM")

for patient in PATIENT_DIST_DICT.keys():
    if patient in PATIENTS:
        print patient
        FILE_OUT1.write("\n"+patient+"\t"+
                        str(PATIENT_DIST_DICT[patient])+"\t"+
                        str(PATIENT_GENE_DICT[patient]) +"\t"+
                        str(PATIENT_ISO_DICT[patient]) +"\t"+
                        str(PATIENT_MUT_DICT[patient]) +"\t"+
                        str(PATIENT_FIRSTISO_DICT[patient]))
FILE_OUT1.close()