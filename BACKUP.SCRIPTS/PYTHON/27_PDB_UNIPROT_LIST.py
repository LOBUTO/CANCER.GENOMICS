#Format pdbtosp.txt to UNICODEB (PDB-UNIPROT ID)

#FORMAT pdbtosp.txt
UNICODE=open("UNIPROT_PDB/pdbtosp.txt")
UNICODEA=UNICODE.readlines()
UNICODE.close()
 
UNICODEA=UNICODEA[25:-5]

X=0
while X!=len(UNICODEA):
    if UNICODEA[X][0]==" ":
        UNICODEA[X]=UNICODEA[X].replace("           ",UNICODEA[X-1][0:11],1)
    X=X+1

UNICODEB=[]
for i in UNICODEA:
    LINE=""
    if len(i)<50 and i[6:11]!="Model":
        LINE=LINE+" "+i[0:4]+" "+i[28:39]+"\n"
    elif len(i)>50 and i[6:11]!="Model":
        LINE=LINE+" "+i[0:4]+" "+i[28:39]+" "+i[50:61]+"\n"
    UNICODEB.append(LINE)

FILE=open("UNIPROT_PDB/UNICODEB","w")
FILE.writelines(UNICODEB)


