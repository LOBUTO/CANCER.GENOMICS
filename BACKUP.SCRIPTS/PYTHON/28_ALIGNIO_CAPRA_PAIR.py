#Make list of EC groups-single member per i.aln (just need the first ID_X or UNIPROT ID per EC groups)

import subprocess 
from Bio import AlignIO

def unique(LIST):
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

#GET i.aln CAPRA FILES THAT HAVE A CORRESPONDING i.sdpL
CAPRA=subprocess.check_output("ls", cwd="capra_bioinf08_data").splitlines()
ALN=filter(lambda a: a[-3:]=="aln", CAPRA) #Use filter from now on
print ALN

SDPA=filter(lambda c: c[-4:]=="sdpA", CAPRA)
print SDPA

ALN_U=[]
for i in ALN:
    if any(i[:-3]==j[:-4] for j in SDPA):
        ALN_U.append(i)
print ALN_U

#GET all files from PDB_LIST_ALL:
PDB_ALL=subprocess.check_output("ls", cwd="PDB_FILES/PDB_LIST_ALL").splitlines()
print PDB_ALL

#GET UNIPROT IDs THAT CONTAIN COMPARABLE PDB:
UNIPROT_UNIQUE=open("UNIPROT_PDB/UNICODEB_UNIPROT").readlines()
UNIPROT=[]
for i in UNIPROT_UNIQUE:
    UNIPROT.append(i[:-1]) #To get rid of "\n"

#CLUSTERING
for i in ALN_U:

    #Join sequences in each i.aln
    ALN1=open("capra_bioinf08_data/%s"%i)
    ALIGNIO=AlignIO.parse(ALN1,"clustal") #Join sequences in each
    
    ALN2=open("BIOPYTHON/ALIGNIO_CAPRA/%sstockholm"%i[:-3],"w") #Save it as i.stockholm
    AlignIO.write(ALIGNIO, ALN2,"stockholm") #Use AlignIO to choose the format, can't read it as it comes out
    ALN1.close() 
    ALN2.close() #Don't forget to close them
    
    ALN3=open("BIOPYTHON/ALIGNIO_CAPRA/%sstockholm"%i[:-3])
    ALN3=ALN3.readlines()
    ALN3=filter(lambda b: b[0]!="#" and b[0]!="/", ALN3) #Get rid of unasable characters # and / in joined sequences file
    
    #Filter out those IDs that do not belong to UNICODEB_UNIPROT OR PDB_LIST_ALL
    ALN4=[]
    for j in ALN3:
        if any (j.split("|")[0][0:4]==k for k in PDB_ALL): #Take if it belongs to the list in the PDB_ALL
            ALN4.append(j)
        elif any(j.split("|")[0]==l for l in UNIPROT): #Take if it belongs to the list in UNIPROT
            ALN4.append(j) 
    
    #Get both EC numbers
    EC=[]
    for m in ALN4:
        EC.append(m.split()[0].split("|")[2])
    EC=unique(EC)
    print EC
    #Accounting for the fact that some may not have more than one element per EC at this point (because not all have a 
    #corresponding PDB file 
    PAR=[]
    if len(EC)>1:        
    
        #Separate into EC lists - At this point the assumption is that there are only two unique numbers in EC
        EC_ZERO=[]
        EC_ONE=[]
        for n in ALN4:
            if n.split()[0].split("|")[2]==EC[0]:
                EC_ZERO.append(n)
            elif n.split()[0].split("|")[2]==EC[1]:
                EC_ONE.append(n)
        
        #Accounting for the fact that some EC groups may contain no elements
        if len(EC_ZERO)>0 and len(EC_ONE)>0:
            PAR.append(EC_ZERO[0]+EC_ONE[0])
            print PAR
            FILE_OUT=open("BIOPYTHON/ALIGNIO_CAPRA_PAIR/%spair"%i[:-3],"w")
            FILE_OUT.writelines(PAR) #WORKS

#413 i.stockholm files that have sdpL files
#137 files will not contain a pair due to lack of PDB file available for comparisson
#276 i.pair files in ALIGNIO_CAPRA_PAIR that contain a single representative of each EC in group 
    