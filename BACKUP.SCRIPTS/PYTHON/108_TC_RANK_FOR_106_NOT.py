#CHECK FOR USAGE OF TC_RANK_V2 FIRST BEFORE USING THIS!!!!!!!!
#WRAPPER FOR TC_RANK SPECIFIC FOR 106.py BASED ON HAVE_SMILES.smi FILES (DATABASES)
#08/17/13
#Dependencies
#    -TC_RANK

#NOT USED IN 106.py ANYMORE!!!!!!
def TC_RANK_FOR_106 (QUERY, SUBJECT, DICT, LOGS): #Takes metabolite NAMES that are present in a personalized
                                                    #smiles dictionary, in the case of 106.py, the dictionary
                                                    #in that script will be based on the HAVE_SMILES.smi
                                                    #WARNING - KEEP IN MIND IF DICT HAS NO KEYS FOR METABOLITES
                                                    #ENTERED THEN TC=0.0
    from FUNCTIONS import TC_RANK
    
    if SUBJECT in DICT and QUERY in DICT:
        #Format files to be taken by TC_RANK
        QUERY_FILE=open(LOGS+"/TRF1.smi", "w")
        QUERY_FILE.write(DICT[QUERY])
        QUERY_FILE.close()
        
        SUBJECT_FILE=open(LOGS+"/TRF2.smi", "w")
        SUBJECT_FILE.write(DICT[SUBJECT]+"\t"+SUBJECT)
        SUBJECT_FILE.close()
        
        #Call TC_RANK
        TC_RANK(LOGS+"/TRF1.smi", LOGS+"/TRF2.smi", LOGS+"/TRF3")
        
        #Call value from file
        FILE1=open(LOGS+"/TRF3")
        FILE_IN1=FILE1.read().splitlines()
        FILE1.close()
        
        TC=float(FILE_IN1[0].split("=")[1])
    
    else:
        TC=0.0 #IF KEY IS NOT FOUND!
    
    return TC

    
#########TESTING#########
import pickle
PICKLE_IN1=open("DATABASES/HMDB_SOURCE/OBJECTS/DICT_METABOLITE_TO_SMILES.pi")
MET_TO_SMILES_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

TEST=TC_RANK_FOR_106("beta-oxycodol",
                "cevimeline N-oxide",
                MET_TO_SMILES_DICT,
                "DATABASES/HMDB_SOURCE/LOGS")
print TEST, type(TEST)
    