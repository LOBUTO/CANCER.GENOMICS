import subprocess
import math
import itertools
from Bio import AlignIO
from Bio.PDB import *

def STOCKHOLM(FILE_ALN, FILE_OUT): #Takes an alignment file of multiple lines per sequence 
                                    #and makes it into a single line per sequence
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
def SPLIT_ALL(OPEN_FILE): #Takes file and returns list of appended splitted() content
    DUMMY=[]
    for i in OPEN_FILE:
        DUMMY.append(i.split())
    return DUMMY
def GET_PDB(PDB_NAME,FOLDER): #PDB FILE - NEEDS Bio.PDB import*, downloads PDB_NAME to specified FOLDER
    pdb1=PDBList()
    pdb1.retrieve_pdb_file(PDB_NAME,pdir=FOLDER)
    return
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
def MIX(LIST1, LIST2,FIX=0): #LIST - MIXES EVERY CONTENT OF LIST1(n) WITH EVERY ELEMENT OF LIST2(m) RETURNING LIST OF mxn ELEMENTS
                            #IF FIX=1 THEN JOINS EVERY ELEMENT OF LIST 2 WITH ALL ELEMENTS OF LIST 1 TO MAKE nxLIST1 LISTS
                            #IF FIX=2 THEN JOINS EVERY ELEMENT OF LIST 1 WITH ALL ELEMENTS OF LIST 2 TO MAKE LIST2xm LISTS
    LIST=[]
    if FIX==0:
        for i in LIST1:
            for j in LIST2:
                LIST.append([i,j])
    elif FIX==1:
        for l in LIST2:
            LIST.append(LIST1+[l])
    elif FIX==2:
        for k in LIST1:
            LIST.append([k]+LIST2)
    
    return LIST

def CLUSTALO_ALIGNIO(FILE_IN, FILE_ALN, FILE_OUT): #ALN FILE and PARSED ALN FILE - Takes in sequences saved to FILE_IN
                                                    # and aligns them using clustal omega to FILE_ALN 
                                                    #and parses them using stockholm, cleans up to FILE_OUT
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
def PDB_RESN_COOR(PDB_FILE, RESN, CHAIN="A"): #LIST - Return atomic coordinate lines from PDB file for residue number within chain of
                                                #interest

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
    
    #Get coordinates that correspond to residue number
    PDBZ=[]
    for f in PDB_IN: #Use chain ATOMS including HETATM embedded if any
        if int(f[22:26])==int(RESN):
            PDBZ.append(f) #list of coordinates of Z atoms
    
    return PDBZ
def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

def PDB_SDPA_LIG(PDB_NAME, SDPA_N, SEQ_FOR_AST, LOGS, CHAIN="A", IONS="TRUE", D=5): #LIST - Returns list of ligands that are within D 
                                                                                    #amstrongs of SDPA_N (Residue #) of interest
                                                                                    #PDB_NAME: name of PDB(4 char), SDPA_N: residue number
                                                                                    #SEQ_FOR_AST:sequence to be used for analysis that
                                                                                    #aligns or has been aligned using MSA to PDB parsed to a
                                                                                    #single line for asterisk replacement method, LOGS: Logs
                                                                                    #folder, CHAIN:chain of interest, IONS: if true will keep
                                                                                    #ions, D:distance threshold
    GET_PDB(PDB_NAME, LOGS)
    
    #Get sequence of PDB file and residue numbering for specific chain of interest
    PDB_AA, PDB_NUM=PDB_SEQ(LOGS+"/pdb%s.ent"%PDB_NAME.lower(), CHAIN)
    
    #Get ligand(s) of PDB
    PDB_LIG_LINES=PDB_LIG(LOGS+"/pdb%s.ent"%PDB_NAME.lower(), IONS)
    
    #Do asterisk replacement of single sequence obtained from multiple sequence alignment (or other source), HAS TO BE SINGLE LINE
    SEQ_AST_IN=SEQ_FOR_AST[0:int(SDPA_N)] + "*" + SEQ_FOR_AST[int(SDPA_N)+1:]
    
    #Save SEQ_AST_IN and PDB_AA to single spaced file
    ALN_IN_FILE=open(LOGS+"/ALN_TEST1","w")
    ALN_IN_FILE.write(PDB_NAME+"_"+CHAIN+" "+SEQ_AST_IN+"\n"+PDB_NAME+" "+PDB_AA)
    ALN_IN_FILE.close()
    
    #Align file and clean up
    CLUSTALO_ALIGNIO(LOGS+"/ALN_TEST1", LOGS+"/ALN_TEST2", LOGS+"/ALN_TEST3")
    
    #Get Z out of ALN_TEST3 file
    FILE_ALN_IN=open(LOGS+"/ALN_TEST3")
    FILE_ALN=SPLIT_ALL(FILE_ALN_IN)
    FILE_ALN_IN.close()
    
    ASTX_SEQ="".join([f for f in itertools.takewhile(lambda x: x!="*", FILE_ALN[0][1])])
    Z=PDB_NUM[len(ASTX_SEQ)]-FILE_ALN[1][1][:len(ASTX_SEQ)].count("-") #Second term is to account for disparities between sequences
    
    #Calculations
    PDB_Z_LINES=PDB_RESN_COOR(LOGS+"/pdb%s.ent"%PDB_NAME.lower(), Z, CHAIN)
    
    PDBL=[] #List of HETATM lines that are withint 5A of Z atoms
    for g in PDB_Z_LINES:
        for q in PDB_LIG_LINES:
            if q[0:6]=="HETATM": #Don't actually need because they are all supposed to be HETATM, just checking
                #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
                DIST=math.sqrt((float(g[30:38])-float(q[30:38]))**2+(float(g[38:46])-float(q[38:46]))**2+(float(g[46:54])-float(q[46:54]))**2)
                if DIST<D and q!=g: #Don't have to put "and q!=g" because I'm already filtering for HETATM
                    PDBL.append("".join(q[17:20].split())) #appending only ligand name    
    
    PDBL=unique(PDBL)
    
    return PDBL
def PDB_RESN(PDB_FILE, RESN, CHAIN="A"): #STR - Returns 3 letter amino acid code corresponding to Residue number (RESN) in chain of interest

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
    
    #Get coordinates that correspond to residue number
    PDB_RES_NAME=[]
    for f in PDB_IN: #Use chain ATOMS including HETATM embedded if any
        if int(f[22:26])==int(RESN):
            PDB_RES_NAME.append(f[17:20]) #list of coordinates of Z atoms
    
    PDB_RES_NAME=unique(PDB_RES_NAME)
    return PDB_RES_NAME[0]
def PDB_SDPA_ZRES(PDB_NAME, SDPA_N, SEQ_FOR_AST, LOGS, CHAIN="A"): #STRING - Returns amino acid three letter code
                                                                                    #corresponding to Z for the SDPA_N given
                                                                                    #PDB_NAME: name of PDB(4 char), SDPA_N: residue number
                                                                                    #SEQ_FOR_AST:sequence to be used for analysis that
                                                                                    #aligns or has been aligned using MSA to PDB parsed to a
                                                                                    #single line for asterisk replacement method, LOGS: Logs
                                                                                    #folder, CHAIN:chain of interest, IONS: if true will keep
                                                                                    #ions, D:distance threshold
    GET_PDB(PDB_NAME, LOGS)
    
    #Get sequence of PDB file and residue numbering for specific chain of interest
    PDB_AA, PDB_NUM=PDB_SEQ(LOGS+"/pdb%s.ent"%PDB_NAME.lower(), CHAIN)
    
    #Do asterisk replacement of single sequence obtained from multiple sequence alignment (or other source), HAS TO BE SINGLE LINE
    SEQ_AST_IN=SEQ_FOR_AST[0:int(SDPA_N)] + "*" + SEQ_FOR_AST[int(SDPA_N)+1:]
    
    #Save SEQ_AST_IN and PDB_AA to single spaced file
    ALN_IN_FILE=open(LOGS+"/ALN_TEST1","w")
    ALN_IN_FILE.write(PDB_NAME+"_"+CHAIN+" "+SEQ_AST_IN+"\n"+PDB_NAME+" "+PDB_AA)
    ALN_IN_FILE.close()
    
    #Align file and clean up
    CLUSTALO_ALIGNIO(LOGS+"/ALN_TEST1", LOGS+"/ALN_TEST2", LOGS+"/ALN_TEST3")
    
    #Get Z out of ALN_TEST3 file
    FILE_ALN_IN=open(LOGS+"/ALN_TEST3")
    FILE_ALN=SPLIT_ALL(FILE_ALN_IN)
    FILE_ALN_IN.close()
    
    ASTX_SEQ="".join([f for f in itertools.takewhile(lambda x: x!="*", FILE_ALN[0][1])])
    Z=PDB_NUM[len(ASTX_SEQ)]-FILE_ALN[1][1][:len(ASTX_SEQ)].count("-") #Second term is to account for disparities between sequences    
    
    #Get amino acid corresponding to Z
    AAA=PDB_RESN(LOGS+"/pdb%s.ent"%PDB_NAME.lower(),Z,CHAIN)
    return AAA

FILE_ENDF=open("LIGAND_OUTPUT/LIG_AA_PAIR_ALT_SINGLE_SET", "w")
    
SEPARATOR=[]

##Get all SDPA files - Got 413
SDPA_IN=subprocess.check_output("ls", cwd="capra_bioinf08_data").splitlines()
SDPA_IN=filter(lambda x: x[-4:]=="sdpA", SDPA_IN) 

#Filter for ones that actually contain SDPAs - Got 240
SDPA_IN=filter(lambda x: len(open("capra_bioinf08_data/%s"%x).readlines())>2, SDPA_IN)

#Work with .aln files
for i in SDPA_IN:
    print i
    #Piece sequences together per .aln file and clean up
    STOCKHOLM("capra_bioinf08_data/%saln"%i[:-4], "LOGS/%ssto"%i[:-4])
    
    #Get pairs of EC numbers per file
    EC_IN=subprocess.Popen(["-c","awk '{print $1}' LOGS/%ssto | cut -d '|' -f3 | sort | uniq"%i[:-4]],
                           stdout=subprocess.PIPE, shell=True) #"-c" tells the shell that the argument is 
                                                                #a command to run, INCLUDE shell=True
    EC=EC_IN.communicate()[0].splitlines()
    
    #Get EC-PDB sets
    STO_IN=open("LOGS/%ssto"%i[:-4])
    STO=SPLIT_ALL(STO_IN)
    STO_IN.close()
    STO=filter(lambda x: len(x[0].split("|")[0])==6, STO)
    
    EC1_PDB, EC2_PDB=[],[] #Store PDB names in each corresponding to EC
    for j in STO:
        if EC[0] in j[0]:
            EC1_PDB.append(j)
        elif EC[1] in j[0]:
            EC2_PDB.append(j)
    
    if len(EC1_PDB)!=0 and len(EC2_PDB)!=0: #Account for sets that are empty (due to UNIPROTs only available in some subgroups)
        #########ADJUST TO ONLY GET ONE REPRESENTATIVE PER SUBGROUP#######
        EC1_PDB=[EC1_PDB[0]]
        EC2_PDB=[EC2_PDB[0]]
        ##################################################################        
        
        PAIR_PDB=MIX(EC1_PDB,EC2_PDB)
        print PAIR_PDB
        #Get SDPAs
        SDPA_IN=open("capra_bioinf08_data/%s"%i)
        SDPA=SPLIT_ALL(SDPA_IN)
        SDPA_IN.close()
        SDPA=SDPA[2:]
        SDPA_N=[]
        for k in SDPA: SDPA_N.extend(k)
        print SDPA_N
        #Work on EC_PDB set of pairs - for each pair in set
        for j in PAIR_PDB:
            FIRST=j[0]
            SECOND=j[1]
            FIRST_LIG_LIST, SECOND_LIG_LIST, FIRST_Z_LIST, SECOND_Z_LIST=[],[],[],[]
            
            #First get a list of all ligands for each one and tag residues with corresponding ligand
            for l in SDPA_N:
                FIRST_LIG=PDB_SDPA_LIG(FIRST[0][0:4], int(l), FIRST[1], "LOGS", CHAIN=FIRST[0][5], IONS="FALSE", D=5)
                if len(FIRST_LIG)!=0:
                    FIRST_LIG_LIST.extend(FIRST_LIG)
                    for F in FIRST_LIG: #Tagging each ligand with each residue
                        FIRST_Z_LIST.append(F+"-"+l) 
                SECOND_LIG=PDB_SDPA_LIG(SECOND[0][0:4], int(l), SECOND[1], "LOGS", CHAIN=SECOND[0][5], IONS="FALSE", D=5)
                if len(SECOND_LIG)!=0:
                    SECOND_LIG_LIST.extend(SECOND_LIG)     
                    for S in SECOND_LIG: #Tagging each ligand with each residue
                        SECOND_Z_LIST.append(S+"-"+l)
            FIRST_LIG_LIST=unique(FIRST_LIG_LIST) 
            SECOND_LIG_LIST=unique(SECOND_LIG_LIST)
            FIRST_Z_LIST=unique(FIRST_Z_LIST)
            SECOND_Z_LIST=unique(SECOND_Z_LIST)
            
            print FIRST_LIG_LIST
            print SECOND_LIG_LIST
            print FIRST_Z_LIST
            print SECOND_Z_LIST
                    
            #Second, iterate over ligands to get intersection of ls
            FIRST_INTS=[]
            for FLL in FIRST_LIG_LIST:
                FIRST_INTS_IN=[FLL]
                for FZL in FIRST_Z_LIST:
                    if FLL in FZL:
                        FIRST_INTS_IN.append(FZL.split("-")[1]) #Make lists of each ligand with its corresponding SDPAs
                FIRST_INTS.append(FIRST_INTS_IN)
            print FIRST_INTS
            
            SECOND_INTS=[]
            for SLL in SECOND_LIG_LIST:
                SECOND_INTS_IN=[SLL]
                for SZL in SECOND_Z_LIST:
                    if SLL in SZL:
                        SECOND_INTS_IN.append(SZL.split("-")[1])
                SECOND_INTS.append(SECOND_INTS_IN)
            print SECOND_INTS
            
            for FI in FIRST_INTS:
                for SI in SECOND_INTS:
                    if len(list(set(FI[1:]) & set(SI[1:])))!=0: #If they both have common SDPAs
                        INT_LIST=list(set(FI[1:]) & set(SI[1:]))
                        print INT_LIST
                        FIRST_AA_LIST, SECOND_AA_LIST=[],[]
                        for IL in INT_LIST:
                            FIRST_AA_LIST.append(PDB_SDPA_ZRES(FIRST[0][0:4],IL, FIRST[1], "LOGS", CHAIN=FIRST[0][5]))
                            SECOND_AA_LIST.append(PDB_SDPA_ZRES(SECOND[0][0:4], IL, SECOND[1], "LOGS", CHAIN=SECOND[0][5]))
                        print FI[0]
                        print "_".join(FIRST_AA_LIST)
                        print SI[0]
                        print "_".join(SECOND_AA_LIST)
                        FILE_ENDF.write(i+"|"+FI[0]+"_"+"_".join(FIRST_AA_LIST)+"|"+SI[0]+"_"+"_".join(SECOND_AA_LIST)+"\n")

FILE_ENDF.close()
            
                        
        
        
        