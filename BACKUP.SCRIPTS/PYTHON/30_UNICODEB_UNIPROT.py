#UNIPROT IDs THAT CONTAINS COMPARABLE PDB
def unique(LIST):
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

PDB_UNIPROT=open("UNIPROT_PDB/UNICODEB").readlines()
UNIPROT=[]
for i in PDB_UNIPROT:
    UNIPROT.extend(i.split()[1:])
UNIPROT=unique(UNIPROT) #Reduced 5-fold to 18840
print len(UNIPROT)
print UNIPROT

FILE_OUT=open("UNIPROT_PDB/UNICODEB_UNIPROT","w")
for i in UNIPROT:
    FILE_OUT.write("%s\n"%i)

FILE_OUT.close()