#HMDB METABOLITES PROCESSING 3
#08/15/13
#Get rid of drugs first

import pickle

#Load dict object
PICKLE_IN1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB.pi")
MET_UNI_DICT=pickle.load(PICKLE_IN1)
PICKLE_IN1.close()

print len(MET_UNI_DICT)
for key in MET_UNI_DICT.iterkeys():
    if "Atorvastatin" in key: print key

#Load drugbank names
DRUG_BANK=open("DATABASES/FDA/FDA_NAMES").read().splitlines() #USING FDA DB INSTEAD OF DRUG BANK NOW!!!!

#Delete keys that have items in DRUG_BANK
COUNT=0
for key_set in MET_UNI_DICT.keys():
    COUNT=COUNT+1
    print COUNT
    if len(set([met.upper() for met in key_set]) & set(DRUG_BANK))>0:
        print key_set
        del MET_UNI_DICT[key_set]

print len(MET_UNI_DICT)

#TEST
for key in MET_UNI_DICT.iterkeys():
    if "Atorvastatin" in key: print key

#Re-store to dict object    
PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB.pi", "w")
pickle.dump(MET_UNI_DICT, PICKLE_OUT1)
PICKLE_OUT1.close()