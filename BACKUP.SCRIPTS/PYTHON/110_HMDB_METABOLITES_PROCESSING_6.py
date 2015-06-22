##HMDB METABOLITES PROCESSING 6 - NOT NEEDED AT THE MOMENT, LOOK AT NOTES!!!! 08/23/13
#08/20/13
#Filtering of proteins in network using the DICT_UNIPROT_TO_SIDE_PRODUCT_PART1.pi and 
#DICT_UNIPROT_TO_SIDE_PRODUCT_PART2.pi dict objects

import pickle

#Load side product dictionaries (from STEP 102.py) {UNIPROT:[SIDE_PRODUCTS]}
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_PART1.pi")
UNI_TO_SP=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

PICKLE_IN2=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_PART2.pi")
UNI_TO_SP.update(pickle.load(PICKLE_IN2))
PICKLE_IN2.close()

print "Number of UNIPROTS IN SIDE PRODUCT DICT", len(UNI_TO_SP)

#Load network dict (from STEP 4.1) {(met1, met2): [UNIPROT1, UNIPROT2]}
PICKLE_IN3=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED.pi")
MET_TO_UNI=pickle.load(PICKLE_IN3)
PICKLE_IN3.close()

#Process
#STEP1: Strip(+/-) on UNI_TO_SP
UNI_TO_SP=dict((X, [Y.strip("+").strip("-") for Y in UNI_TO_SP[X]]) for X in UNI_TO_SP)    

#STEP2: Get total number of uniprots in MET_TO_UNI dict
MET_TO_UNI_UNIPROTS=[item for subset in MET_TO_UNI.values() for item in subset]

ALL_LIGANDS=set([item.upper() for subset in MET_TO_UNI.keys() for item in subset])
SD_LIGANDS=set([item.upper() for subset in UNI_TO_SP.values() for item in subset])
print "Number of ligands in SD",len(SD_LIGANDS)
print "Number of all ligands", len(ALL_LIGANDS)
print len(SD_LIGANDS)
print len(SD_LIGANDS & ALL_LIGANDS)

#STEP3: Filter UNI_TO_SP with it:
PRE_TOTAL=len(set(MET_TO_UNI_UNIPROTS))
Z=0
UNI_COUNT=0
for key in MET_TO_UNI.keys():
    PRE=len(MET_TO_UNI[key]) #To see difference
    
    MET_TO_UNI[key]=[X for X in MET_TO_UNI[key] if UNI_TO_SP.has_key(X)==False 
                     or len(set([met.upper() for met in key]) & set([sp.upper() for sp in UNI_TO_SP[X]]))==0]
    
    POST=len(MET_TO_UNI[key]) #To see difference
    
    if POST<PRE:#To give me account of how much it is being optimized
        Z=Z+1
        UNI_COUNT=UNI_COUNT+(PRE-POST)
        print key

POST_TOTAL= len(set([item for subset in MET_TO_UNI.itervalues() for item in subset]))

print "PRE_TOTAL", PRE_TOTAL
print "POST_TOTAL", POST_TOTAL
print "AFFECTED ENTRIES", Z
print "DELETED NUMBER OF UNIPROT NODES", UNI_COUNT

PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED_PLUS.pi", "w")
pickle.dump(MET_TO_UNI, PICKLE_OUT1)
PICKLE_OUT1.close()