#FORMAT pdbtosp.txt FOR PDB->UNIPROT OBJECT (RECENT)
import pickle

FILE_IN=open("DATABASES/pdbtosp_201308.txt").read().splitlines()

#First replace blank pdb identifiers
X=1
while X!=len(FILE_IN):
    if FILE_IN[X][0]==" ":
        FILE_IN[X]=FILE_IN[X].replace("           ", FILE_IN[X-1][0:11],1)

    X=X+1

#Then get pdb and uniprot only and get rid of MODELS
PDB_UNIPROT_LINES=[]
for record in FILE_IN:
    if len(record)<50 and record[6:11]!="Model":
        PDB_UNIPROT_LINES.append(record[0:4]+" "+record[41:47])
    elif len(record)>50 and record[6:11]!="Model":
        PDB_UNIPROT_LINES.append(record[0:4]+" "+record[41:47]+" "+record[63:69])

#Make dictionary # SAVE DICT AS OBJECT WHERE KEY IS PDB AND VALUES ARE LISTS OF AVALIBALE UNIPROTS
PDB_DICT=dict((x[0:4],[]) for x in PDB_UNIPROT_LINES)
for record in PDB_UNIPROT_LINES:
    PDB_DICT[record.split()[0].strip()]=PDB_DICT[record.split()[0].strip()]+record.split()[1:]

PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_PDB_TO_UNIPROT_ALL.pi", "w")
pickle.dump(PDB_DICT, PICKLE_OUT1)