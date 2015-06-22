#01/20/14
#QUICK AND DIRTY PROCESS THOSE CANCERS THAT HAVE SYNONYMOUS MUTATIONS
#SEMINORMALIZE BY GENE LENGTH

import subprocess, collections,pickle

MISSENSE_FILES=filter(lambda x: "mut" in x, subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/by_cancer_aa_info").splitlines())
MISSENSE_FILES=["DATABASES/CANCER_DATA/by_cancer_aa_info/"+X for X in MISSENSE_FILES]

SILENT_FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/synonymousMuts").splitlines()
SILENT_FILES=["DATABASES/CANCER_DATA/synonymousMuts/"+X for X in SILENT_FILES]

#Load GENE_LENGTH dictionary
PICKLE_IN1=open("DATABASES/UNIPROT/OBJECTS/012014_UNIPROT_GENE_LENGTH.pi")
DICT_GENE_LEN=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

def PROCESS_MUT (FILES):
    
    OUT_DICT=dict((X.split("/")[-1].split("_")[0],{}) for X in FILES) #To get name of CANCER, each CANCER is a key
    PATIENT_DICT=dict((X.split("/")[-1].split("_")[0],0) for X in FILES)
    
    #Load files
    for record in FILES:
        FILE_IN=open(record)
        RECORD=[Y.split("\t") for Y in FILE_IN]
        FILE_IN.close()
        
        #Get TOTAL NUMBER OF PATIENTS
        TOTAL_PATIENTS=set([C[1] for C in RECORD])
        
        #Make gene:#isoform count dict 
        ISOFORM_COUNT=list(set([Y[0] for Y in RECORD]))
        ISOFORM_COUNT=dict(collections.Counter([Z.split(".")[0] for Z in ISOFORM_COUNT]))

        #Make dict of total mutations per gene (gene:mutationsxpatients)
        GENE_COUNT=dict((Z,0) for Z in ISOFORM_COUNT.keys())
        
        for flat in RECORD:
            GENE_COUNT[flat[0].split(".")[0]]=GENE_COUNT[flat[0].split(".")[0]] + 1
        
        #Store to dict as cancer:gene:mutaionsxpatients - Normalized to number of isoforms and SEMINORMALIZED BY GENE LENGTH
        OUT_DICT[record.split("/")[-1].split("_")[0]]=dict((W, 
                                                            float(GENE_COUNT[W])/(ISOFORM_COUNT[W])) 
                                                           for W in GENE_COUNT.keys())
        
        #Store to dict patient SET
        PATIENT_DICT[record.split("/")[-1].split("_")[0]]=TOTAL_PATIENTS
        
    return (OUT_DICT, PATIENT_DICT)

#Call function
MISSENSE_DICT, MISSENSE_PATIENT_DICT=PROCESS_MUT(MISSENSE_FILES)
SILENT_DICT, SILENT_PATIENT_DICT=PROCESS_MUT(SILENT_FILES)
CANCERS=SILENT_DICT.keys() #Either one would do

#Write out to file
FILE_OUT1=open("DATABASES/CANCER_DATA/LOG/012114_CANCER_PROCESSED_COUNTS", "w")
FILE_OUT1.write("CANCER"+"\t"+"GENE"+"\t"+"SILENT_COUNT"+"\t"+"MISSENSE_COUNT"+"\t"+"REST")

for cancer in CANCERS:
    CANCER_GENES=list(set(SILENT_DICT[cancer].keys() + MISSENSE_DICT[cancer].keys()))
    ALL_CANCER_PATIENTS= len(MISSENSE_PATIENT_DICT[cancer] | SILENT_PATIENT_DICT[cancer])
    
    for gene in CANCER_GENES:
        
        if DICT_GENE_LEN.has_key(gene): #NEED TO HAVE LENGTH TO CALCULATE "REST"
            
            FILE_OUT1.write("\n"+cancer+"\t"+gene.upper()) #KEEP IN MIND THE UPPER()
            
            #First write SILENT count for gene
            if SILENT_DICT[cancer].has_key(gene):
                SILENT=SILENT_DICT[cancer][gene]
                FILE_OUT1.write("\t"+ str(SILENT))
            else:
                SILENT=1
                FILE_OUT1.write("\t"+str(SILENT)) #USING 1 FOR COMPARISSON!!!
            
            #Then MISSENSE count
            if MISSENSE_DICT[cancer].has_key(gene):
                MISSENSE=MISSENSE_DICT[cancer][gene]
                FILE_OUT1.write("\t"+ str(MISSENSE))
            else:
                MISSENSE=1
                FILE_OUT1.write("\t"+str(MISSENSE)) #USING 1 FOR COMPARISSON!!!!!
            
            #Then REST count
            REST=ALL_CANCER_PATIENTS*DICT_GENE_LEN[gene] - (SILENT+MISSENSE)
            FILE_OUT1.write("\t"+str(REST))
            
FILE_OUT1.close()
        


        
        
        
        
            