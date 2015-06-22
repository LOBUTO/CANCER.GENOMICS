#HMDB METABOLITES PROCESSING 7B
#08/22/13
#Making file for Cytoscape
#Have to be tab_delimited files

import pickle, csv

#Load file
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB_FILTERED.pi")
METABO_TO_UNIPROTS=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

#Get rid of water
for key_set in METABO_TO_UNIPROTS.keys():
    if "Water" in key_set:
        del METABO_TO_UNIPROTS[key_set]

#Write to 2 Tab separated file, one for interaction per line and another for attribute
CSV_OUT1=open("NETWORK/TSF_FILES/HMDB_FILTERED.csv", "w")
CSV_OUT2=open("NETWORK/TSF_FILES/HMDB_FILTERED_ATR.csv", "w") #Get node attributes

for key_set in METABO_TO_UNIPROTS.keys():
    CSV_OUT2.write("_".join(min(key_set, key=len).split())+"\t"+"METABOLITE"+"\n")
    
    for UNIPROT in METABO_TO_UNIPROTS[key_set]:
        CSV_OUT1.write("_".join(min(key_set, key=len).split())+"\t"+UNIPROT+"\n")
        CSV_OUT2.write(UNIPROT+"\t"+"UNIPROT"+"\n")

CSV_OUT1.close()
CSV_OUT2.close()