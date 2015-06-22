#MAPPING FUNCTION FROM PDB_RES TO UNIPROT_RES - USE UNIPROT SITE TO TEST
#6/29/13 - Relies on EMBOSS asterisk replacement method
#Dependencies:
#    GET_PDB
#    PDB_SEQ
#    UNIPROT_SEQ
#    EMBOSS_ALIGNIO
#    PDB_CHAIN_TO_UNIPROT_CLASS
#INTRODUCE AS DICT FORMAT: {PDB1|CHAIN:[RES1, RES2,...],....} => PDB_DICT
#RESULTS ARE IN THE DICT FORMAT: {PDB1|CHAIN=UNIPROT: [[URES1, NAM], [URES2, NAM],...], PDB2|CHAIN=UNIPROT: [.....  }
#CAUTION - PDB|X that does not have a UNIPROT match will return to DICT in the format {PDB|X:"None"}!!!!

#FIX!!!! CANNOT USE CLUSTALO ANYMORE, USE EMBOSS NEEDLE INSTEAD    
def PDB_RES_TO_UNIPROT_RES (PDB_DICT, LOG_FOLDER): #THROUGH *REPLACEMENT METHOD
    from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, EMBOSS_ALIGNIO, PDB_CHAIN_TO_UNIPROT_CLASS
    import pickle
    import urllib2
    import itertools
    
    #Import Biolip-based object that maps pdb-chain to uniprot #FORMAT IS {pdb|chain:uniprot, ......}
    #DICT_FILE1=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_CHAIN_TO_UNIPROT.pi?attredirects=0&d=1")
    #DICT_PDB_CHAIN_TO_UNIPROT=pickle.load(DICT_FILE1)

    UNIPROT_DICT={}
    print PDB_DICT
    for record in PDB_DICT.iterkeys():
        print record
        #Parse each record
        PDB=record.split("|")[0].upper()
        CHAIN=record.split("|")[1].upper()
        RES_LIST=PDB_DICT[record]
        
        #Get pdb file - Will download as pdb****.ent where **** is the 4 digit pdb identifier
        PDB_FILE=GET_PDB(PDB.upper(), LOG_FOLDER)
        
        #Get sequence and numbering from file
        PDB_AA, PDB_NUM=PDB_SEQ(LOG_FOLDER+"/pdb%s.ent"%PDB.lower(), CHAIN.upper())
        
        #Get Uniprot using PDB_CHAIN_TO_UNIPROT_CLASS
        UNIPROT=PDB_CHAIN_TO_UNIPROT_CLASS([PDB.upper()+"|"+CHAIN.upper()]) #Takes LIST of PDB|X, in this case just one in list
        WITH_UNIPROT=UNIPROT.with_uniprots() 
        WITHOUT_UNIPROT=UNIPROT.without_uniprots()#Those without found uniprot will return {PDB|X:"None"}    
        
        #Check if uniprot is found and continue
        if len(WITH_UNIPROT)>0:
            
            #Get Uniprot sequence
            UNIPROT_AA=UNIPROT_SEQ(WITH_UNIPROT[record.upper()])
            
            #Do asterisk replacement method - THIS WILL DO A EMBOSS PER RESIDUE, IF EVER THERE IS A BETTER WAY TO DO THIS, IMPLEMENT HERE
            UNIPROT_RES=[]
            for res in RES_LIST:
                if res in PDB_NUM: #IMPORTANT - To account for the fact that there are false pdb numbering in CSA tables                
                    #Replace with asterisk at PDB residue number
                    PDB_AST_IND=PDB_NUM.index(res)
                    PDB_AST_AA=PDB_AA[:PDB_AST_IND]+"*"+PDB_AA[PDB_AST_IND+1:]
                    
                    #Align pdb and uniprot sequences
                    PDB_FASTA=open(LOG_FOLDER+"/PRTUR1.fasta", "w") #Write sequences PER file for EMBOSS processing
                    UNIPROT_FASTA=open(LOG_FOLDER+"/PRTUR2.fasta", "w")
                    PDB_FASTA.write(">"+PDB+"\n"+PDB_AST_AA)
                    UNIPROT_FASTA.write(">"+WITH_UNIPROT[record.upper()]+"\n"+UNIPROT_AA)
                    PDB_FASTA.close()
                    UNIPROT_FASTA.close()    
                    EMBOSS_ALIGNIO(LOG_FOLDER+"/PRTUR1.fasta", LOG_FOLDER+"/PRTUR2.fasta", LOG_FOLDER+"/LOG3", LOG_FOLDER)
                    
                    #Get uniprot residue number from file
                    FILE_IN1=open(LOG_FOLDER+"/LOG3")
                    ALN_SEQUENCES=FILE_IN1.read().splitlines()
                    PDB_AL=ALN_SEQUENCES[0].split()[1]
                    UNI_AL=ALN_SEQUENCES[1].split()[1]
                    FILE_IN1.close()
                    
                    AST_COUNT=len([x for x in itertools.takewhile(lambda x: x!="*", PDB_AL)])
                    UNI_FIX=UNI_AL[:AST_COUNT].count("-") #In case of missmaches in uniprot sequence
                    UNI_COUNT=AST_COUNT+1-UNI_FIX
                    
                    UNIPROT_RES.append([UNI_COUNT,UNI_AL[AST_COUNT]])
        
            UNIPROT_DICT[record.upper()+"="+WITH_UNIPROT[record.upper()]]=UNIPROT_RES
        
        #If no uniprot found for PDB|X, will update as {PDB|X:"None"}
        else:
            UNIPROT_DICT.update(WITHOUT_UNIPROT)
        
    return UNIPROT_DICT 

SAMPLE={"9PAP|A":[175,19,25,159], "1bzy|c":[133, 134, 169, 137,104]}

TEST=PDB_RES_TO_UNIPROT_RES(SAMPLE, "LOGS")
print "TEST", TEST