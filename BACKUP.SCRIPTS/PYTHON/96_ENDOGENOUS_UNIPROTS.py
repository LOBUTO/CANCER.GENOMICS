#GET ALL ENDOGENOUS UNIPROTS FROM NETWORK
from FUNCTIONS import unique

FILE_IN=open("NETWORK/ENDOGENOUS.ADJ").read().splitlines()

ENDOGENOUS_UNIPROT=[]
for record in FILE_IN:
    if len(record.split())>1:
        ENDOGENOUS_UNIPROT=ENDOGENOUS_UNIPROT+[X.split("|")[1].strip() for X in record.split()[1:]]
print len(ENDOGENOUS_UNIPROT), ENDOGENOUS_UNIPROT
ENDOGENOUS_UNIPROT=unique(ENDOGENOUS_UNIPROT)
print len(ENDOGENOUS_UNIPROT)


FILE_OUT=open("NETWORK/ENDOGENOUS_UNIPROT", "w")
for record in ENDOGENOUS_UNIPROT:
    FILE_OUT.write(record+"\n")
FILE_OUT.close()