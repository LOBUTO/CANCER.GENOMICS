#1/6/13
#Take hypergeometric results and filter out patients that have significantly enriched metabolic mutated genes
#FILES USED 
#    "NETWORK/FROM_R_ANALYSIS/010614_SIGNIFICANT_METABOLITES"
#    "NETWORK/PMVs/010714.pi"
#PRODUCED: Tables per cancer for heatmap analysis in R
import pickle

#Load table
FILE_IN1=open("NETWORK/FROM_R_ANALYSIS/012014_SIGNIFICANT_METABOLITES")
SIG_CANCER_MET=[X.split("\t") for X in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

#Extract significant metabolites per cancer - Create dictionary cancer:{sig_mets}
CANCERS=dict((X[1],[]) for X in SIG_CANCER_MET)

for record in SIG_CANCER_MET:
    CANCERS[record[1]]=CANCERS[record[1]]+[record[0]]

#Load PMVs - cancer:patient:metabolite:#mutated genes
PICKLE_IN1=open("DATABASES/CANCER_DATA/OBJECTS/121713_DICT_PMV_NETWORK100.pi")
PMV=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

for cancer in PMV.keys():
    if cancer not in ["LUSC","READ", "GBM","KIRC", "UCEC", "OV", "BRCA","COAD"]:
        del PMV[cancer]

#Filter out patients and metabolites that do not posses significance
for cancer in PMV.keys():
    print cancer
    
    for patient in PMV[cancer].keys():
        
        #If patients do not have significant metabolites remove them
        if len(set(PMV[cancer][patient].keys()) & set(CANCERS[cancer]))==0:
            del PMV[cancer][patient]
        
        #Otherwise only keep metabolites that are affected per patient
        else:
            for metabolite in PMV[cancer][patient].keys():
                if metabolite not in CANCERS[cancer]:
                    del PMV[cancer][patient][metabolite]


#Write out table for R_ANALYSIS per CANCER (So that each table is comprised of only metabolites that affect a particular cancer)
for cancer in PMV.keys():
    print "TABLE_%s"%cancer
    
    FILE_OUT1=open("NETWORK/R_ANALYSIS/012014_SIGNIFICANT_%s_PMV_TABLE"%cancer, "w")

    FILE_OUT1.write("PATIENT")
    
    #Write a column per metabolite in each cancer
    for met in CANCERS[cancer]:
        FILE_OUT1.write("\t"+met)
        
    #Write number of affected genes per metabolite per patient
    for patient in PMV[cancer].keys():
        
        FILE_OUT1.write("\n"+patient)
        
        for met in CANCERS[cancer]: #Granted is writing in the same order as column order
            
            #If patient does not have cancer significant metabolite write zero
            if met not in PMV[cancer][patient].keys():
                FILE_OUT1.write("\t"+"0")
            
            #Otherwise write number of mutated proteins associated with metabolite
            else:
                FILE_OUT1.write("\t"+ str(PMV[cancer][patient][met]))
    
    FILE_OUT1.close()
                    
#PROCEED TO HEATMAP AND PCA IN R  
        
        