#HMDB METABOLITES PROCESSING 4.1
#08/19/13
#Filter and process synonyms and TC_RANK filtered sets

import pickle
from FUNCTIONS import TC_RANK_V2

#Load dict with keys as tuple of names and values as lists of uniprot_ids
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB.pi")
MET_UNI_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

#Load MET_SMILES_DICT for TC_RANK_V2 function
PICKLE_IN2=open("DATABASES/HMDB_SOURCE/OBJECTS/DICT_METABOLITE_TO_SMILES.pi")
MET_TO_SMILES_DICT=pickle.load(PICKLE_IN2)
PICKLE_IN2.close()

#Process - Keep in mind keys are tuples that have to be transiently converted to sets for data manipulation
FILTERED_MET_UNI_DICT={}
COUNT=None
while COUNT!=len(FILTERED_MET_UNI_DICT): #TO ITERATE OVER UNTIL FILTERED DICT HAS NO MORE SYNONYMS OR TC MATCHES AT ALL!
    COUNT=len(FILTERED_MET_UNI_DICT)
    for key_set in MET_UNI_DICT.keys(): #keys() instead of iterkeys() creates a copy that can call for "del"
                                        #of keys when necessary in the dictionary, BUT KEEP IN MIND THAT
                                        #EVEN IF A KEY THAT HAS BEEN NOT ITERATED FOR IS DELETED, IT WILL
                                        #STILL ITERATE OVER IT AND WILL RETURN AN ERROR IF THEY VALUE OF THE
                                        #KEY IS CALLED FOR
        
        if key_set in MET_UNI_DICT: #To account for the fact that iteration will call key_set even though "key"
                                    #has potentially been removed from dictionary
            
            print key_set
            
            #Store all that have same key in key_sets in dictionary, including itself as list of sets. Intersection
            MATCHING_KEY_SETS=[record for record in MET_UNI_DICT.keys() if len(set(key_set) & set(record))>0]
            #Store as list of tuples(as they were) instead of list of sets.
            #The reasoning for this is that even if it doesn't have another record, it will find itself and
            #store only one record. We need to keep it as a tuple so that the key is recognized in the order
            #the elements are in it so it does get removed. This record will then go on through the pipeline
            #through frozenset (which takes tuples, even if the list contains only one of them).
            #print MATCHING_KEY_SETS 
            
            #Store all records that have a record[0] of TC=1.0 when compared to key_set[0]
            #This step compares the first element of key_set(name of metabolite) against all first left_over metabolites
            #of MET_UNI_DICT ([0]). Stores those key_sets (record) that have their name in list of names of matching TC
            TC_RANKING=TC_RANK_V2(key_set[0], [Y[0] for Y in MET_UNI_DICT.keys()], MET_TO_SMILES_DICT, "DATABASES/HMDB_SOURCE/LOGS",
                                                                                             FILTER=1.0)    
            MATCHING_TC=[record for record in MET_UNI_DICT.keys() if record[0] in TC_RANKING] #Store if name matches scored metabolites
            
            #UPDATE TO MATCHING_KEY_SETS
            #This is assuming that if AxB=1.0 and BxC=1.0 that AxC. Because if we are iterating over a deleted
            #entry then that would mean that MATCHING_TC would also be empty along with MATCHING_KEY_SETS
            MATCHING_KEY_SETS=MATCHING_KEY_SETS+MATCHING_TC
            MATCHING_KEY_SETS=list(set(X for X in MATCHING_KEY_SETS)) #Make whole list of sets unique
            
            #Process
            if len(MATCHING_KEY_SETS)!=0:
            
                #Combine all into single key
                NEW_KEY_SET=tuple(frozenset().union(*MATCHING_KEY_SETS)) # asterisks(*) breaks internal sets (or tuples)
                
                #Combine all uniprot ids as list of sets - keep in mind actual key_sets are tuples
                ALL_UNIPROT_LISTS=[set(MET_UNI_DICT[KEYS]) for KEYS in MATCHING_KEY_SETS]
                ALL_UNIPROT=list(frozenset().union(*ALL_UNIPROT_LISTS))
                
                #Update to new dictionary. We could use current one because of two reasons:
                #1. If it the same single record. If we remove if before we add it, it won't go to the end of the line
                    #in the dict iteration because we are iterating over a copy of original call. So it will be 
                    #perfectly stored
                #2. If it is a combination of records, the old records will be removed, which is fine, and
                    #the new record will be added without being iterated over later on, again, because we are
                    #iterating over a copy of the original call
                    
                FILTERED_MET_UNI_DICT[NEW_KEY_SET]=ALL_UNIPROT
                
                #Remove all used entries - don't need to do it outside loop since empty MATCHING_KEY_SETS
                #would mean that the record does not exist anymore and does not need to be removed.
                for used_key in MATCHING_KEY_SETS:
                    del MET_UNI_DICT[used_key]
        
        #Give a count of how much is left
        print len(MET_UNI_DICT), len(FILTERED_MET_UNI_DICT)
        FILE_OUT=open("DATABASES/HMDB_SOURCE/LOGS/TRAININGV2", "w")
        FILE_OUT.write(str(len(MET_UNI_DICT))+" "+ str(len(FILTERED_MET_UNI_DICT))+"\n")
        FILE_OUT.close()
    
print len(FILTERED_MET_UNI_DICT)

#Store as new dict object

PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED.pi", "w")
pickle.dump(FILTERED_MET_UNI_DICT, PICKLE_OUT1)
PICKLE_OUT1.close()