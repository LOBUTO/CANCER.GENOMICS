#SCRIPT 3 - Introducing PFAM delimited sequences to IDs in clusters

import subprocess
import urllib,urllib2
from Bio.PDB import *
import itertools

def UNIPROT_SEQ(UNIPROT): #STRING - Feed UNIPROT ID and returns FASTA sequence
    SEQ=""
    url=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta"%UNIPROT).readlines()
    for i in url[1:]:
        SEQ=SEQ+i[:-1]
    return SEQ

def MAKE_DICT_LIST(FILE_IN): #DICTIONARY - Makes GO Term dictionary with OUT_FILE_DIC
    FILE=open(FILE_IN)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split())
    
    print FILEA
    
    MOL_DICT={} #VERY FAST - DICTIONARY
    for key,val in FILEA:
        MOL_DICT.setdefault(key,[]).append(val)
    
    FILE.close()
    return MOL_DICT

def GET_PDB(PDB_NAME,FOLDER): #PDB FILE - NEEDS Bio.PDB import*, downloads PDB_NAME to specified FOLDER
    pdb1=PDBList()
    pdb1.retrieve_pdb_file(PDB_NAME,pdir=FOLDER)
    return

#Assuming these are the files we have

def PDB_SEQ(PDB_FILE, CHAIN="A"): #STRING + LIST(int) - Give PDB file and chain and returns sequence as string and residue 
                                    #numbering as a list of integers
    #AMINO ACID DICTIONARY
    AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     #modified amino acids
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'ASG', 
     'CME':'C', 'CSX':'C','CSO':'C','CSS':'C', 'SEP':'S','TPO':'T','LLP':'K',
     'SUI':'X', "NEP":"H","FTR":"W","TRQ":"W","TQQ":"W", "TTQ":"W","1TQ":"W", 
     "0AF":"W", "TRW":"W", "TOQ":"W" }
        
    PDB_IN=open(PDB_FILE).readlines()
    
    #Get rid of comments
    PDB_IN= PDB_IN[len([a for a in itertools.takewhile(lambda x: x[0:4]!="ATOM",PDB_IN)]):]
    
    #Account for unmentioned chain "_"
    if CHAIN=="_":
        CHAIN=PDB_IN[0][21:22]
    
    #Get rid of ANISOU lines and get ATOMS for chain of interest including HETATM within chain
    PDB_IN=filter(lambda x:x[:6]!="ANISOU" and x[21:22]==CHAIN,PDB_IN)
    PDB_IN=[b for b in itertools.takewhile(lambda y: y[0:3]!="TER",PDB_IN)]
    
    #Get resname and resnumber sequences
    RES=[]
    RES_NUM=[]
    RES_NAM=""
    for i in PDB_IN:
        if i[17:20]+i[22:26] not in RES: 
            RES.append(i[17:20]+i[22:26])
            RES_NUM.append(int(i[22:26]))
            RES_NAM=RES_NAM+AADICT[i[17:20]]
    
    print PDB_IN[0] #To check
    print PDB_IN[-1] #To check
    return RES_NAM, RES_NUM

SEP=[]

# STEP 1 - Get ID_PFAM-LIST for mini dictionaries
CLUSTERS=subprocess.check_output("ls", cwd="CLASS_PROJECT/GROUPS_FILTER").splitlines()

GO_TAGS=[]
for i in CLUSTERS:
    TAG=subprocess.Popen(["less", "CLASS_PROJECT/GROUPS_FILTER/%s"%i], stdout=subprocess.PIPE)
    GO_TAGS.extend(TAG.communicate()[0].splitlines())
    
ID_FOR_CLUSTER=open("CLASS_PROJECT/ID_FOR_CLUSTER","w") #Found 509 unique
for i in GO_TAGS:
    ID_FOR_CLUSTER.write(i.split()[1]+" "+i.split()[2]+"\n")
ID_FOR_CLUSTER.close()

#STEP 2 - Add PFAM domains
#Create parsable PDB_PFAM1 and SWISS_PFAM lists
subprocess.call("awk '{print $1  $2, $3, $4}' CLASS_PROJECT/PDB_PFAM1 |sort > temp1",shell=True)
subprocess.call("awk '{print $1  $2, $3, $4}' CLASS_PROJECT/SWISS_PFAM |sort > temp2",shell=True)

#Join by ID_PFAM
for i in CLUSTERS:
    subprocess.call("sort -k2.1 CLASS_PROJECT/GROUPS_FILTER/%s | awk '{print $1, $2  $3}' | sort -k2.1 > temp"%i, shell=True)
    subprocess.call("join -o '1.1 1.2 2.2 2.3' -1 2 temp temp1 > temp3", shell=True)
    subprocess.call("join -o '1.1 1.2 2.2 2.3' -1 2 temp temp2 >> temp3", shell=True)
    subprocess.call("sort temp3 > CLASS_PROJECT/GROUPS_FILTER_PFAM/%s"%i, shell=True)

#Clean up files
ID_PFAM=subprocess.check_output("ls", cwd="CLASS_PROJECT/GROUPS_FILTER_PFAM").splitlines()
for i in ID_PFAM:
    U_FILE=[]
    U_FILE_IN=open("CLASS_PROJECT/GROUPS_FILTER_PFAM/%s"%i)
    for j in U_FILE_IN:
        U_FILE.append(j.split())
    U_FILE_IN.close()
    C_FILE_OUT=open("CLASS_PROJECT/GROUPS_T/%s"%i, "w")
    for k in U_FILE:
        C_FILE_OUT.write(k[0]+"|"+k[1][0:(len(k[1])-7)]+"|"+k[2]+"-"+k[3]+"\n")
    C_FILE_OUT.close()

#STEP 3- Get sequences:
#Get SWISS-UNIPROT minidictionary
SWISS_UNIPROT=MAKE_DICT_LIST("CLASS_PROJECT/SWISS_UNIPROT_FOR_CLUSTER")

#Call for UNIPROT SEQUENCE
ID_PFAM=subprocess.check_output("ls", cwd="CLASS_PROJECT/GROUPS_FILTER_PFAM").splitlines()
for i in ID_PFAM:
    CLUSTER_FILE=[]
    CLUSTER_FILE_IN=open("CLASS_PROJECT/GROUPS_FILTER_PFAM/%s"%i)
    for j in CLUSTER_FILE_IN:
        CLUSTER_FILE.append(j.split("|"))
    CLUSTER_FILE_IN.close()
    print CLUSTER_FILE
    
    CLUSTER_FILE_OUT=open("CLASS_PROJECT/GROUPS_FILTER_SEQ/%s"%i,"w")
    for k in CLUSTER_FILE:
        if SWISS_UNIPROT.has_key(k[1]):
            UNIPROT=SWISS_UNIPROT[k[1]][0]
            print UNIPROT
            CLUSTER_FILE_OUT.write(k[0]+"|"+k[1]+"|"+k[2][:-1]+" "+UNIPROT_SEQ(UNIPROT)[int(k[2].split("-")[0]):int(k[2].split("-")[1])]+"\n")
        
        else:
            GET_PDB(k[1][0:4],"CLASS_PROJECT/LOGS")
            PDB_AA, PDB_NUM=PDB_SEQ("CLASS_PROJECT/LOGS/pdb%s.ent"%k[1][0:4].lower(), CHAIN=k[1][5])
            PFAM_LEFT=len([n for n in itertools.takewhile(lambda y: y!=int(k[2].split("-")[0]), PDB_NUM)])
            PFAM_RIGHT=len([m for m in itertools.takewhile(lambda z: z!=int(k[2].split("-")[1]), PDB_NUM)])
            CLUSTER_FILE_OUT.write(k[0]+"|"+k[1]+"|"+k[2][:-1]+" "+PDB_AA[PFAM_LEFT:PFAM_RIGHT]+"\n")
    CLUSTER_FILE_OUT.close()


