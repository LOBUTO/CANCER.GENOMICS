#HMDB METABOLITES PROCESSING 7A
#08/20/13
#Finalize processing of METABO_TO_UNIPROT dict object into a readable NETWORKX file
#NO MORE BETA STEP NEEDED AT THE MOMENT - LOOK AT NOTES 8/23/13!!!

import pickle
"""
#######TEMPORAL STEP "BETA" - ELIMINATE ALL COFACTORS/METAL&IONS######
#Load side product dictionaries (from STEP 102.py) {UNIPROT:[SIDE_PRODUCTS]} and get rid off (+/-)
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_ALL_PART1.pi")
UNI_TO_SP=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

for FN in range(2,10):
    PICKLE_IN_FN=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_ALL_PART%s.pi"%str(FN))
    UNI_TO_SP.update(pickle.load(PICKLE_IN_FN))
    PICKLE_IN_FN.close()

#PICKLE_IN2=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_PART2.pi")
#UNI_TO_SP.update(pickle.load(PICKLE_IN2))
#PICKLE_IN2.close()


UNI_TO_SP=dict((X, [Y.strip("+").strip("-") for Y in UNI_TO_SP[X]]) for X in UNI_TO_SP)
##################
"""
#Load network dict (from STEP 4.1) {(met1, met2): [UNIPROT1, UNIPROT2]}
PICKLE_IN3=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED.pi")
METABO_TO_UNIPROTS=pickle.load(PICKLE_IN3)
PICKLE_IN3.close()
"""
#################
#Eliminate all cofactors from Network
SP_NAMES=set([SP.upper() for SP_SET in UNI_TO_SP.values() for SP in SP_SET])
for record in METABO_TO_UNIPROTS.keys():
    if len(set([X.upper() for X in record]) & SP_NAMES)>0:
        del METABO_TO_UNIPROTS[record]

print len(METABO_TO_UNIPROTS)
#################
"""
#Write to file as is: (Make to encode utf-8 for NETWORKX reading and join metabolite name by "_")
FILE_OUT1=open("NETWORK/MEATBOLITE_TO_UNIPROT_HMDB_FILTERED.ADJ", "w")
for MET_SET, UNI_LIST  in METABO_TO_UNIPROTS.items():
    FILE_OUT1.write(u"%s"%"_".join(min(MET_SET, key=len).split()).encode("utf-8")) #Shortest name of met_set chosen
    for uniprot in UNI_LIST:
        FILE_OUT1.write(" "+u"%s"%uniprot.encode("utf-8"))
    FILE_OUT1.write("\n")
FILE_OUT1.close()