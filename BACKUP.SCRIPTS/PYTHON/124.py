#NETWORK MODELING #1

#010713 - Obtain KEGG pathways for KEGG ids
#PRODUCED - "KEGG_FILES/010714_DICT_PATH_CO.pi" - KEGG_PATH:["DESCRIPTION"],["COMPOUNDS"]

from urllib2 import urlopen
import pickle

#Obtain human pathways in KEGG
URL_IN1=urlopen("http://rest.kegg.jp/list/pathway/hsa").read().splitlines()
KEGG_HUMAN_PATH=[X.split("\t") for X in URL_IN1]

#Build PATH:COMPOUNDS dict
KEGG_DICT=dict((X[0], {}) for X in KEGG_HUMAN_PATH)

for record in KEGG_HUMAN_PATH:
    KEGG_DICT[record[0]]["DESCRIPTION"]=record[1]
    
for pathway in KEGG_DICT.keys():
    
    MAP="map"+pathway.split(":")[1][3:]
    URL_IN2=urlopen("http://rest.kegg.jp/link/cpd/" + MAP).read().splitlines()
    
    try:
        COMPOUNDS=[Y.split("\t")[1].split(":")[1] for Y in URL_IN2]
    
        KEGG_DICT[pathway]["COMPOUNDS"]=COMPOUNDS

    except:
        print "Exception",MAP
        #Remove those that have no COMPOUNDS
        del KEGG_DICT[pathway]
        pass 


#Pickle out
PICKLE_OUT1=open("KEGG_FILES/010714_DICT_PATH_CO.pi", "w")
pickle.dump(KEGG_DICT, PICKLE_OUT1)
PICKLE_OUT1.close()
