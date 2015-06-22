#SCRIPT 6 - Validate SDPA
from Bio.PDB import *
import itertools
import subprocess
from Bio import AlignIO
import math

def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

def GET_PDB(PDB_NAME,FOLDER): #PDB FILE - NEEDS from Bio.PDB import*, downloads PDB_NAME to specified FOLDER
    pdb1=PDBList()
    pdb1.retrieve_pdb_file(PDB_NAME,pdir=FOLDER)
    return

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
     'SUI':'X', "NEP":"H", "FTR":"W", "TRQ":"W", "TQQ":"W", "TTQ":"W","1TQ":"W", 
     "0AF":"W", "TRW":"W" ,"TOQ":"W" }
        
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

def PDB_LIG(PDB_FILE, IONS="TRUE"): #LIST - Takes PDB file and returns lines corresponding to LIGANDS, TRUE argument holds if
                                    #ions are to be included
    PDB_IN=open(PDB_FILE).readlines()
    
    #Filter for HETATMs after last chain
    PDB_IN=[b for b in itertools.takewhile(lambda y:y[0:3]!="TER",PDB_IN[::-1])]
    PDB_IN=PDB_IN[::-1]
    
    #Filter for ANISOU and CONNECT
    PDB_IN=filter(lambda w: w[0:6]=="HETATM",PDB_IN)
    
    #Filter for HOHs:
    PDB_IN=filter(lambda z:z[17:20]!="HOH",PDB_IN)
    
    #ION conditional
    if IONS=="TRUE":
        return PDB_IN
    elif IONS=="FALSE":
        PDB_IN=filter(lambda a: len("".join(a[17:20].split()))>2,PDB_IN)
        return PDB_IN

def PDB_LIG_TRUE(PDB_FILE, IONS="TRUE"): #TRUE or FALSE - Give PDB and returns whether or not PDB has ligands, IONS=TRUE argument
                                            #holds unless ions do not want to be included
    PDB_IN=open(PDB_FILE).readlines()
    
    #Filter for HETATMs after last chain
    PDB_IN=[b for b in itertools.takewhile(lambda y:y[0:3]!="TER",PDB_IN[::-1])]
    PDB_IN=PDB_IN[::-1]
    
    #Filter for ANISOU and CONNECT
    PDB_IN=filter(lambda w: w[0:6]=="HETATM",PDB_IN)
    
    #Filter for HOHs:
    PDB_IN=filter(lambda z:z[17:20]!="HOH",PDB_IN)
    
    #ION conditional
    if IONS=="TRUE":
        PDB_IN=PDB_IN
    elif IONS=="FALSE":
        PDB_IN=filter(lambda a: len("".join(a[17:20].split()))>2,PDB_IN)
    
    #Testing
    if len(PDB_IN)>1: return "TRUE"
    else: return "FALSE"    

def PDB_COOR(PDB_FILE, CHAIN="A"): #LIST - Return atomic coordinate lines from PDB file for chain of interest

    #AMINO ACID DICTIONARY
    AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     #modified amino acids
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'ASG', 
     'CME':'C', 'CSX':'C','CSO':'C','CSS':'C', 'SEP':'S','TPO':'T','LLP':'K',
     'SUI':'X', "NEP":"H", "FTR":"W", "TRQ":"W", "TQQ":"W", "TTQ":"W","1TQ":"W", 
     "0AF":"W", "TRW":"W" ,"TOQ":"W" }
        
    PDB_IN=open(PDB_FILE).readlines()
    
    #Get rid of comments
    PDB_IN= PDB_IN[len([a for a in itertools.takewhile(lambda x: x[0:4]!="ATOM",PDB_IN)]):]
    
    #Account for unmentioned chain "_"
    if CHAIN=="_":
        CHAIN=PDB_IN[0][21:22]
    
    #Get rid of ANISOU lines and get ATOMS for chain of interest including HETATM within chain
    PDB_IN=filter(lambda x:x[:6]!="ANISOU" and x[21:22]==CHAIN,PDB_IN)
    PDB_IN=[b for b in itertools.takewhile(lambda y: y[0:3]!="TER",PDB_IN)]
    return PDB_IN

def CLUSTALO_ALIGNIO(FILE_IN, FILE_ALN, FILE_OUT): #ALN FILE and PARSED ALN FILE - Takes in sequences and aligns them using
                                                    #clustal omega, returns aligned file (FILE_ALN) and parsed file (FILE_OUT)
                                                    
    subprocess.check_output(["CLUSTALO/clustalo", "-i", FILE_IN, "-o", FILE_ALN, "--outfmt=clu", "--force"])
    
    FILE=open(FILE_ALN)
    alignments=AlignIO.parse(FILE, "clustal") #Parses the aligned file to join sequences
    OUTPUT=open(FILE_OUT, "w")
    AlignIO.write(alignments, OUTPUT, "stockholm") #Saves the parsed file in the stockholm format
    
    FILE.close()
    OUTPUT.close()

    #Get rid of "#" and "/" containing lines
    PARSED_IN=open(FILE_OUT)
    PARSED=PARSED_IN.readlines()
    PARSED_IN.close()
    PARSED=filter(lambda x: x[0]!="#" and x[0]!="/", PARSED)
    PARSED_OUT=open(FILE_OUT, "w")
    PARSED_OUT.writelines(PARSED)
    PARSED_OUT.close()

def PDB_RES_LIG(PDB, SEQ, LOGS, CHAIN="A",D=5, IONS="TRUE"): #LIST - Returns numbers mapping to fed sequence of residues
                                                                #that are within D amstrongs of ligands
                                                                #PDB is ID of interest, SEQ is sequence obtained from alignment
                                                                #parsed through CLUSTALO_ALIGNIO, LOGS is the LOG folder, D is
                                                                #restrictive distance, IONS=TRUE would filter out all ligands of
                                                                #length < 3 
    GET_PDB(PDB, LOGS) #WARNING - SET YOUR PREFERRED PDB FOLDER IF YOU ARE INTERESTED IN KEEPING PDB
    RES_AA, RES_N=PDB_SEQ(LOGS+"/pdb%s.ent"%PDB.lower(), CHAIN)
    LIG=PDB_LIG(LOGS+"/pdb%s.ent"%PDB.lower(), IONS)
    
    #Map SEQ TO FILE:
    TEMP1=open(LOGS+"/temp1", "w") #WARNING - WILL SAVE ALL "UNSEEN" FILES TO LOGS FOLDER
    TEMP1.write(PDB+"_"+CHAIN+" "+SEQ+"\n"+PDB+" "+RES_AA)
    TEMP1.close()
    
    CLUSTALO_ALIGNIO(LOGS+"/temp1", LOGS+"/temp2", LOGS+"/temp3") #WARNING - WILL SAVE ALL "UNSEEN" FILES TO LOGS FOLDER    
    
    #Get coordinate lines corresponding to SEQ
    ALN=open(LOGS+"/temp3")
    ALN_A=[]
    for i in ALN: ALN_A.append(i.split())
    print ALN_A
    LEFT=len([j for j in itertools.takewhile(lambda x: x=="-", ALN_A[0][1])])
    RIGHT=len(ALN_A[0][1]) - len([l for l in itertools.takewhile(lambda y: y=="-", ALN_A[0][1][::-1])])
    PDB_LEFT=RES_N[LEFT]
    PDB_RIGHT=RES_N[RIGHT]
    
    PDB_LINES=PDB_COOR(LOGS+"/pdb%s.ent"%PDB.lower(),CHAIN)
    PDB_LINES=filter(lambda z: int(z[22:26])>=PDB_LEFT and int(z[22:26])<PDB_RIGHT, PDB_LINES)
    
    #Get residues in PDB within distance constraint (D)
    Z=[]
    for o in PDB_LINES:
        for p in LIG:
            DIST=math.sqrt((float(o[30:38])-float(p[30:38]))**2+
                           (float(o[38:46])-float(p[38:46]))**2+
                           (float(o[46:54])-float(p[46:54]))**2)
            if DIST<D:
                Z.append(o[22:26])
    Z=unique(Z)
    print Z
    
    #Map back to SEQ
    Z_SEQ=[]
    for q in Z:
        Z_SEQ.append(ALN_A[0][1][RES_N.index(int(q))])
    print Z_SEQ
    
    Z_SEQ_N=[]
    for r in Z:
        UC=RES_N.index(int(r))-LEFT
        C=UC+SEQ[0:UC].count("-")
        Z_SEQ_N.append(C)
    return Z_SEQ_N
    
def PDB_SDPA_TEST(PDB, SEQ, SDPA_FILE, LOGS, CHAIN="A", D=5, IONS="TRUE"):#LIST - returns list of SDPA-Y/N depending on whether
                                                                            #it passes the D constraint
    GET_PDB(PDB, LOGS) #WARNING - SET YOUR PREFERRED PDB FOLDER IF YOU ARE INTERESTED IN KEEPING PDB
    RES_AA, RES_N=PDB_SEQ(LOGS+"/pdb%s.ent"%PDB.lower(), CHAIN)
    LIG=PDB_LIG(LOGS+"/pdb%s.ent"%PDB.lower(), IONS)    
    
    #Get SDPAs
    SDPA_IN=open(SDPA_FILE)
    SDPA=[]
    for i in SDPA_IN:
        SDPA.extend(i.split())

    #Calculate per SDPA:
    SDPA_RESULTS=[]
    for j in SDPA:
        SEQ_AST=SEQ[0:int(j)]+"*"+SEQ[int(j)+1:]
        
        #Map SEQ TO FILE:
        TEMP1=open(LOGS+"/temp1", "w") #WARNING - WILL SAVE ALL "UNSEEN" FILES TO LOGS FOLDER
        TEMP1.write(PDB+"_"+CHAIN+" "+SEQ_AST+"\n"+PDB+" "+RES_AA)
        TEMP1.close()
        
        #Align
        CLUSTALO_ALIGNIO(LOGS+"/temp1", LOGS+"/temp2", LOGS+"/temp3") #WARNING - WILL SAVE ALL "UNSEEN" FILES TO LOGS FOLDER       
        
        #Get Z
        ALG_IN=open(LOGS+"/temp3")
        ALG=[]
        for l in ALG_IN: ALG.append(l.split())
        ALG_IN.close()
        
        Z=RES_N[ALG[0][1].index("*")]
        PDB_LINES=PDB_COOR(LOGS+"/pdb%s.ent"%PDB.lower(), CHAIN)
        PDB_LINES=filter(lambda x: int(x[22:26])==Z, PDB_LINES)
        
        RESIDUE=[]
        for m in PDB_LINES:
            for n in LIG:
                DIST=math.sqrt((float(m[30:38])-float(n[30:38]))**2+
                               (float(m[38:46])-float(n[38:46]))**2+
                               (float(m[46:54])-float(n[46:54]))**2)
                if DIST<D:
                    RESIDUE.append(m[22:26]+"-"+m[17:20]+" "+n[17:20])
        RESIDUE=unique(RESIDUE)
        print RESIDUE
        if len(RESIDUE)>0:
            SDPA_RESULTS.append(j+"-"+"Y")
        else:
            SDPA_RESULTS.append(j+"-"+"N")
    
    return SDPA_RESULTS

SEPARATOR=[]

#FIRST - Find files that contain SDPAs out of 17 clusters analyzed - FOUND 8 Files containig SDPAs
CLUSTER_IN=subprocess.Popen(["find", "CLASS_PROJECT/SDPA", "-mindepth", "1", "!", "-size","0"], stdout=subprocess.PIPE)
CLUSTER=CLUSTER_IN.communicate()[0].splitlines()
print CLUSTER

#SECOND - Test parsed files - Obtained 237 entries
SDPA_TEST=open("CLASS_PROJECT/RESULTS/SDPA_TEST", "w")
for i in CLUSTER:
    PARSED_IN=open("CLASS_PROJECT/ALIGNMENT_PARSED/%s"%i.split("/")[2][:-5])
    PARSED=[]
    for j in PARSED_IN: PARSED.append(j.split()) 
    PARSED_IN.close()
    
    #Test each PDB in file (Can't test UNIPROT since their association to specific PDBs is ambigous)
    for k in PARSED:
        if len(k[0].split("|")[1])==6:
            SDPA_TEST.write(k[0].split("|")[1]+"|"+k[0].split("|")[2]+" ")
            TEST=PDB_SDPA_TEST(k[0].split("|")[1][0:4], k[1], i ,"CLASS_PROJECT/LOGS",
                                CHAIN=k[0].split("|")[1][5], D=5, IONS="FALSE")
            for l in TEST:
                SDPA_TEST.write(l+" ")
            SDPA_TEST.write("\n")

SDPA_TEST.close()

#THIRD - Get Gold standard - Obtained 237 entries
GOLD_SD=open("CLASS_PROJECT/RESULTS/GOLD_SD","w")
for i in CLUSTER:
    PARSED_IN=open("CLASS_PROJECT/ALIGNMENT_PARSED/%s"%i.split("/")[2][:-5])
    PARSED=[]
    for j in PARSED_IN: PARSED.append(j.split()) 
    PARSED_IN.close()
    
    #Get TRUE residues from each PDB in file (Can't test UNIPROT since their association to specific PDBs is ambigous)
    for k in PARSED:
        if len(k[0].split("|")[1])==6:
            GOLD_SD.write(k[0].split("|")[1]+"|"+k[0].split("|")[2]+" ")
            GOLD=PDB_RES_LIG(k[0].split("|")[1][0:4],k[1],"CLASS_PROJECT/LOGS", 
                             CHAIN=k[0].split("|")[1][5], D=5,IONS="FALSE")
            for l in GOLD:
                GOLD_SD.write(str(l)+"-Y"+" ")
            GOLD_SD.write("\n")

GOLD_SD.close()
            
            
