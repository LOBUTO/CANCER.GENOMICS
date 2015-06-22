#1/6/14
#Produce table of metabolites for hypergeometric calculation in R
#Network used: HMDB_NETWORK_TC100FILTER.pi 
#Produced table for R_ANALYSIS

import pickle

#Load network DICT[MET]["GENES"]
PICKLE_IN1=open("DATABASES/HMDB_SOURCE/OBJECTS/NETWORKS_TC_FILTERS/HMDB_NETWORK_TC100FILTER.pi")
NETWORK100=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

NETWORK100_GENES=list(set([X.strip().upper() for Y in NETWORK100.keys() for X in NETWORK100[Y]["GENES"]]))

#Load patient dictionary of the form {cancer:patient:gene:#mutations} 
PICKLE_IN2=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PATIENT_MUTATED_GENES.pi")
PATIENT_DICT=pickle.load(PICKLE_IN2)
PICKLE_IN2.close()

#CREATE TABLE
FILE_OUT1=open("NETWORK/R_ANALYSIS/121813_METABOLITE_TABLE_MUTATIONS", "w")

FILE_OUT1.write("METABOLITE" +"\t"+ "CANCER"+"\t"+"TOTAL_MET_MUTATIONS"+"\t"+"TOTAL_NON_MET_MUTATIONS"+
                "\t"+"CANCER_MET_MUTATIONS"+"\t"+"CANCER_NON_MET_MUTATIONS")

count=1 #to output count 
for metabolite in NETWORK100.keys():
    print count
    
    METABOLITE_GENES=NETWORK100[metabolite]["GENES"]
    
    #Get number of mutations across all cancers related to metabolite (WHITE BALLS)
    METABOLITE_MUTATIONS={}
    for cancer in PATIENT_DICT.keys():
        METABOLITE_MUTATIONS[cancer]=[mutations for Y in PATIENT_DICT[cancer].keys() for gene,mutations in PATIENT_DICT[cancer][Y].items() 
                              if gene in METABOLITE_GENES] 
        METABOLITE_MUTATIONS[cancer]=int(round(sum(METABOLITE_MUTATIONS[cancer]))) #round because hypergeometric calculation only takes integers
    
    #Get number of mutations across all cancer that are not related to metabolite(BLACK BALLS)
    NON_METABOLITE_MUTATIONS={}
    for cancer in PATIENT_DICT.keys():
        NON_METABOLITE_MUTATIONS[cancer]=[mutations for Y in PATIENT_DICT[cancer].keys() for gene,mutations in PATIENT_DICT[cancer][Y].items() 
                              if gene not in METABOLITE_GENES] 
        NON_METABOLITE_MUTATIONS[cancer]=int(round(sum(NON_METABOLITE_MUTATIONS[cancer]))) #round because hypergeometric calculation only takes integers
        
    #Write to table
    for cancer in PATIENT_DICT.keys():
        
        FILE_OUT1.write("\n" + metabolite +"\t" +
                        cancer + "\t" +
                        str(sum(METABOLITE_MUTATIONS.values())) + "\t" +
                        str(sum(NON_METABOLITE_MUTATIONS.values())) + "\t" +
                        str(int(round(METABOLITE_MUTATIONS[cancer]))) + "\t" +
                        str(int(round(NON_METABOLITE_MUTATIONS[cancer])))
                        )
    
    count=count+1

FILE_OUT1.close()
            