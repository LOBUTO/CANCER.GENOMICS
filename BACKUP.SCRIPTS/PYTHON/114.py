#HMDB METABOLITES PROCESSING 8
#08/22/13
#Convert Bi-partite network into unipartite protein network using jaccard and store as
#Cytoscape readable file for A TAB B TAB WEIGHT

import pickle

#Load file - Again, unfiltered with all cofactors network
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED.pi")
METABO_TO_UNIPROTS=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

#Get rid of water
for key_set in METABO_TO_UNIPROTS.keys():
    if "Water" in key_set:
        del METABO_TO_UNIPROTS[key_set]

#Get all uniprots
ALL_UNIPROTS=list(set([uni for UNIPROT_LIST in METABO_TO_UNIPROTS.values() for uni in UNIPROT_LIST]))
print len(ALL_UNIPROTS)

#Send to UNIPROT.org for service