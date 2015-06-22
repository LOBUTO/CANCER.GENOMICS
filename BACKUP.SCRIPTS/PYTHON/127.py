    #01/20/14
#Process Pawel's 1000G data
#Getting total count of mutations per gene (mutations x patients)
#Non-synonymous vs synonymous per patient (Normalized to number of isoforms)
#SEMINORMALIZED BY GENE_LENGTH -NOT YET

import subprocess, collections,pickle

FOLDER1=subprocess.check_output("ls", cwd="DATABASES/1000G_MUT").splitlines()

FILE_OUT1=open("DATABASES/CANCER_DATA/1000GENOME/012014_PROCESSED_COUNTS", "w")
FILE_OUT1.write("GENE"+ "\t"+ "SILENT"+"\t"+"MISSENSE"+"\t"+"REST")

#Load gene_length dictionary
PICKLE_IN1=open("DATABASES/UNIPROT/OBJECTS/012014_UNIPROT_GENE_LENGTH.pi")
DICT_GENE_LEN=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

TOTAL_PATIENTS=2184.0

for record in FOLDER1:
    FILE_IN=open("DATABASES/1000G_MUT/%s"%record)
    THOUSAND=[X.split("\t") for X in FILE_IN.read().splitlines()]
    THOUSAND_GENES=list(set([GENE[0].split(".")[0] for GENE in THOUSAND])) #All genes in project
    FILE_IN.close()
    
    #Divide into SILENT(synonymous) and MISSENSE (non-synonymous)
    SILENT=filter(lambda x: x[2]!="Missense_Mutation", THOUSAND)
    MISSENSE=filter(lambda x: x[2]!="Silent", THOUSAND)
    
    #Get dictionary of number of isoforms per gene for SILENT and MISSENSE - GENE:# of isoforms
    ISOFORM_DICT_SILENT=list(set([X[0] for X in SILENT]))
    ISOFORM_DICT_SILENT=dict(collections.Counter([Y.split(".")[0] for Y in ISOFORM_DICT_SILENT]))
    
    ISOFORM_DICT_MISSENSE=list(set([X[0] for X in MISSENSE]))
    ISOFORM_DICT_MISSENSE=dict(collections.Counter([Y.split(".")[0] for Y in ISOFORM_DICT_MISSENSE]))
    
    #Create patient count dict - Number of mutations per gene - GENE: #Patientsxmutations
    PATIENT_DICT_SILENT=dict((X,0) for X in     )
    PATIENT_DICT_MISSENSE=dict((X,0) for X in ISOFORM_DICT_MISSENSE)
    
    for line in THOUSAND:
        if line[2]=="Silent":
        
            PATIENT_DICT_SILENT[line[0].split(".")[0]]=PATIENT_DICT_SILENT[line[0].split(".")[0]]+int(line[1].split("-")[1])
            
        else:
            PATIENT_DICT_MISSENSE[line[0].split(".")[0]]=PATIENT_DICT_MISSENSE[line[0].split(".")[0]]+int(line[1].split("-")[1])
    
    #Write to file - normalize by number of isoforms and SEMINORMALIZE BY PROTEIN LENGTH
    for gene in THOUSAND_GENES:
    
        if DICT_GENE_LEN.has_key(gene.upper()):    
            FILE_OUT1.write("\n"+gene.upper()) #KEEP IN MIND UPPER()!!
            
            #SILENT count first
            if ISOFORM_DICT_SILENT.has_key(gene):
                SILENT=float(PATIENT_DICT_SILENT[gene]) / (ISOFORM_DICT_SILENT[gene] )
                FILE_OUT1.write("\t" + str(SILENT) )
            else:
                SILENT=1
                FILE_OUT1.write("\t" + str(SILENT)) #TYPING 1 FOR COMPARISSON!!!!!
            
            #Then MISSENSE
            if ISOFORM_DICT_MISSENSE.has_key(gene):
                MISSENSE=float(PATIENT_DICT_MISSENSE[gene]) / (ISOFORM_DICT_MISSENSE[gene]) 
                FILE_OUT1.write("\t" + str(MISSENSE) ) 
            else:
                MISSENSE=1
                FILE_OUT1.write("\t" + str(MISSENSE)) #TYPING 1 FOR COMPARISSON!!!!
                
            #Then REST
            REST=DICT_GENE_LEN[gene.upper()]*TOTAL_PATIENTS-(SILENT+MISSENSE)
            FILE_OUT1.write("\t"+str(REST))
            
        else:
            print record, gene
        
FILE_OUT1.close()
    
            