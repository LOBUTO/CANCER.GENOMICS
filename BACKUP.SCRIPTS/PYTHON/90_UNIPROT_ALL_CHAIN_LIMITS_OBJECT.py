from Bio import SwissProt
import pickle

HANDLE=open("DATABASES/uniprot_sprot.dat")

DICT={}
for record in SwissProt.parse(HANDLE):
    CHAIN_VALUES=[]
    AC=[]
    for ac in record.accessions:
        AC.append(ac) 
    for feature in record.features:
        if feature[0]=="CHAIN":
            CHAIN_VALUES=CHAIN_VALUES+[str(feature[1]), str(feature[2])]
    print CHAIN_VALUES
    
    if any(X.isdigit()==False for X in CHAIN_VALUES) or not CHAIN_VALUES:
        CHAIN_START=1
        CHAIN_END=record.sequence_length
    else:
        CHAIN_START=min(int(X) for X in CHAIN_VALUES)
        CHAIN_END=max(int(X) for X in CHAIN_VALUES)
    print AC[0]
    print CHAIN_START, CHAIN_END
    DICT[AC[0]]=[int(CHAIN_START), int(CHAIN_END)]

PICKLE_OUT=open("DATABASES/OBJECTS/DICT_UNIPROT_CHAIN_LIMITS.pi", "w")
pickle.dump(DICT,PICKLE_OUT)