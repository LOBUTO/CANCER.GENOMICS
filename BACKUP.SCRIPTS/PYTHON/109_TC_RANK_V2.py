#TC_RANK IMPLEMENTATION
#08/19/13
#Dependencies:
#    Have Babel installed properly, make sure path is correct
#Returns list of molecules that pass filter

def TC_RANK_V2(QUERY, LIST, DICT, LOGS, FILTER=1.0, FP="FP2"): #QUERY = Name of metabolite to be analyzed
                                                #LIST = List of metabolites to be analyzed against                                        
                                                #DICT = Dictionary of MET_TO_SMILES for all entries
                                                #LOGS = Logs
                                                #FILTER= higher than it will be displayed
                                                #Returns LIST of molecules names that pass filter, or empty list if none,
                                                #if no key found in dictionary provided, it will also return an 
                                                #empty list
    import subprocess
    
    #Check first if query is in dictionary
    if DICT.has_key(QUERY):
    
        #Prep files
        QUERY_FILE=open(LOGS+"/TCRV21.smi", "w")
        QUERY_FILE.write(DICT[QUERY])
        QUERY_FILE.close()
        
        LIST_FILE=open(LOGS+"/TCRV22.smi", "w")
        for record in LIST:
            if DICT.has_key(record): #Checking for keys (metabolites or key_sets) that do not actually have smiles code
                LIST_FILE.write(DICT[record]+"\t"+record+"\n")
        LIST_FILE.close()
        
        #CHECK THAT LIST_FILE  HAS AT LEAST ONE RECORD IN IT TO RUN WITH, OTHERWISE RETURN EMPTY LIST
        CHECK_LIST=subprocess.check_output(["wc", "-l", "%s/TCRV22.smi"%LOGS])
        COUNT=int(CHECK_LIST.split()[0].strip())
        
        if COUNT>0:
        
            #Run Babel
            TERMINAL_RESULT=subprocess.Popen(["-c", "/usr/local/bin/babel %s/TCRV21.smi %s/TCRV22.smi -ofpt -xf%s"%(LOGS, LOGS, FP)], 
                                             stdout=subprocess.PIPE, shell=True)
            TERMINAL_RESULT=TERMINAL_RESULT.communicate()[0].splitlines()
        
            #Filter out troublesome lines
            TERMINAL_RESULT=filter(lambda x:len(x)>3 and x!="Possible superstructure of first mol", TERMINAL_RESULT)
            
            #Rank in descending order
            RANK=[[float(X.split("=")[1].strip()), X.split("Tanimoto")[0].strip().strip(">")] for X in TERMINAL_RESULT]
            
            RANK=sorted(RANK, reverse=True)
            
            #Return molecule names that pass filter
            FILTER=[Y[1] for Y in RANK if Y[0]>=1.0]
        
        else: #IF LIST_FILE HAS NO COMPATIBLE KEYS AT ALL IT WILL RETURN AN EMPTY LIST!!!
            FILTER=[]
    
    else: #IF QUERY IS NOT IN DICTIONARY PROVIDED IT WILL RETURN AN EMPTY LIST!!!
        FILTER=[]
        
    return FILTER


#####TESTING#######
import pickle

PICKLE_IN1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB.pi")
MET_UNI_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

PICKLE_IN2=open("DATABASES/HMDB_SOURCE/OBJECTS/DICT_METABOLITE_TO_SMILES.pi")
MET_TO_SMILES_DICT=pickle.load(PICKLE_IN2)
PICKLE_IN2.close()

#Check that dictionary is not empty, (empty only if there are no smiles whatsoever for molecules)
if len(MET_TO_SMILES_DICT)!=0:
    #Test one triglyceride against all of the endogenous unfiltered
    CALL=TC_RANK_V2(('TG(18:1(11Z)/14:1(9Z)/20:1(11Z))', '1-(11Z-Octadecenoyl)-2-(9Z-tetradecenoyl)-3-(11-eicosenoyl)-glycerol', '1-Vaccenoyl-2-myristoleoyl-3-eicosenoyl-glycerol', 'TAG(18:1/14:1/20:1)', 'TAG(52:3)', 'TG(18:1/14:1/20:1)', 'TG(52:3)', 'Tracylglycerol(18:1/14:1/20:1)', 'Tracylglycerol(52:3)', 'Triacylglycerol', 'Triglyceride')[0]
                , [X[0] for X in MET_UNI_DICT.keys()], MET_TO_SMILES_DICT, "DATABASES/HMDB_SOURCE/LOGS")

RESULT=[X for X in MET_UNI_DICT.keys() if X[0] in CALL]

print len(CALL)
print len(RESULT)