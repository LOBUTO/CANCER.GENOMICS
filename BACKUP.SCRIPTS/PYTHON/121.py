#12/18/13
#Use patient dictionary {cancer:patient:mutated genes:#mutations} to construct metabolic vectors
#Network used: HMDB_NETWORK_TC100FILTER.pi 
#Produced PMV of the form {cancer:patient:metabolite:normalized count}

import pickle, collections

#Load patient dictionary of the form {cancer:patient:gene:#mutations} 
PICKLE_IN1=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PATIENT_MUTATED_GENES.pi")
PATIENT_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

#Load network DICT[MET]["GENES"]
PICKLE_IN2=open("DATABASES/HMDB_SOURCE/OBJECTS/NETWORKS_TC_FILTERS/HMDB_NETWORK_TC100FILTER.pi")
NETWORK100=pickle.load(PICKLE_IN2)
PICKLE_IN2.close()

NETWORK100_GENES=list(set([X.strip().upper() for Y in NETWORK100.keys() for X in NETWORK100[Y]["GENES"]]))

#Create output dict
DICT=dict((cancer, {}) for cancer in PATIENT_DICT.keys())

#Construct vectors showing normalized counts of how many proteins interact with metabolites
for cancer in PATIENT_DICT.keys():
    
    CANCER_PATIENTS=PATIENT_DICT[cancer].keys()
    CANCER_GENES=set([gene for patient in PATIENT_DICT[cancer].keys() for gene in PATIENT_DICT[cancer][patient].keys()])
    CANCER_METABOLITES=set([X for X in NETWORK100.keys() if len(set(NETWORK100[X]["GENES"]) & CANCER_GENES)>0 ])
    
    #####INFO#####
    print cancer, len(CANCER_PATIENTS), len(CANCER_GENES), len(CANCER_GENES & set(NETWORK100_GENES)), len(CANCER_METABOLITES)
    ##############
    
    for patient in CANCER_PATIENTS:
        
        PATIENT_GENES=[X.strip().upper() for X in PATIENT_DICT[cancer][patient].keys()]
        PATIENT_METABOLIC_GENES=list(set(PATIENT_GENES) & set(NETWORK100_GENES))
        
        if len(PATIENT_METABOLIC_GENES)>0: #ACCOUNT FOR THE FACT THAT SOME PATIENT MAY NOT HAVE ANY MUTATED METABOLIC GENES
            
            #Get all metabolite items per patients counted by number of mutated metabolic genes
            PATIENT_METABOLITES=[]
            for GENE in PATIENT_METABOLIC_GENES:
                for MET in NETWORK100.keys():
                    if GENE in NETWORK100[MET]["GENES"]:
                        PATIENT_METABOLITES=PATIENT_METABOLITES+[MET]
            
            #Dict of MET:#MUTATED_GENES
            PATIENT_METABOLITES_DICT=dict(collections.Counter(PATIENT_METABOLITES))
            
            #Normalize dict to total number of metabolites per patient
            PATIENT_METABOLITES=dict((MET, float(NUMBER)/sum(PATIENT_METABOLITES_DICT.values()) ) for MET,NUMBER in PATIENT_METABOLITES_DICT.items())
            
            #Store to main DICT
            DICT[cancer][patient]=PATIENT_METABOLITES

#Pickle_out of the form {cancer:patient:metabolite:normalized count}
PICKLE_OUT1=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PMV_NETWORK100.pi", "w")
pickle.dump(DICT, PICKLE_OUT1)
PICKLE_OUT1.close()

for i,j in DICT["KIRP"]["TCGA-HE-7129-01A"].items():
    print i,j
