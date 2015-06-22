#12/17/13
#Convert processed TCGA files into patient objects
#Object produced {cancer:patient:gene:#mutations}

import subprocess, pickle

#Get files
TCGA_FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS").splitlines()

#Create output dict
DICT=dict((cancer.split(".")[0], {}) for cancer in TCGA_FILES)

#Break into patient vectors
for record in TCGA_FILES:
    print record
    
    FILE_IN1=open("DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/"+record)
    CANCER_RECORD=[X.split("\t") for X in FILE_IN1.read().splitlines()]
    FILE_IN1.close()
    
    CANCER_PATIENTS=list(set([Y[1] for Y in CANCER_RECORD]))
    
    #List patient:{mutated genes:#mutations}
    PATIENT_DICT=dict((patient,{}) for patient in CANCER_PATIENTS)
    
    for patient in CANCER_PATIENTS:
        
        PATIENT_RECORD=filter(lambda REC: REC[1]==patient, CANCER_RECORD)
        PATIENT_ISOFORMS=list(set([ISO[0] for ISO in PATIENT_RECORD])) #unique isoforms
        PATIENT_GENES=[GENE.split(".")[0] for GENE in PATIENT_ISOFORMS] #with duplicates (each count per gene is a mutation)
        
        #We define the number of mutations per gene as the number of mutations over all isoforms of a gene divided by the number of isoform for that gene
        PATIENT_DICT[patient]=dict((gene, float(PATIENT_GENES.count(gene))/ [ISO_GENE.split(".")[0] for ISO_GENE in PATIENT_ISOFORMS].count(gene)) 
                                   for gene in list(set(PATIENT_GENES)))
     
    DICT[record.split(".")[0]]=PATIENT_DICT

#Pickle out final form {cancer:patient:gene:#mutations} 
PICKLE_OUT1=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PATIENT_MUTATED_GENES.pi", "w")
pickle.dump(DICT, PICKLE_OUT1)
PICKLE_OUT1.close()

print DICT["KIRP"]["TCGA-HE-7129-01A"]
        