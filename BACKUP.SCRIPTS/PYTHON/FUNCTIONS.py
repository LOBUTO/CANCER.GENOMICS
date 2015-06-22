import os
import math
import subprocess
import itertools
import urllib,urllib2
from Bio.PDB import *
from Bio import AlignIO

def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

def FILT_LIG(LIST): #LIST - FILTERS FOR LIGANDS IN LIST
    LIST=filter(lambda x:len(x)<4,LIST)
    return LIST

def FILT_RES(LIST): #LIST- FILTER FOR ANYTHING BESIDES LIGANDS IN LIST
    LIST=filter(lambda x:len(x)>4,LIST)
    return LIST    

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
 
def SMI_FILE(MOL,FILE_IN): #SAVES SMILE FILE, MOL IS MOLECULE OF INTEREST AND FILE_IN IS FILE TO STORE IT IN
    SMILES=open("Components-smiles-oe.smi").readlines() #UPDATE THIS IF SMILES DATABASE IS PLACED SOMEWHERE ELSE
    for i in SMILES:
        if i.split()[1]==MOL:
            FILE=open(FILE_IN,"w")
            FILE.writelines(" ".join(i.split()))
            FILE.close()

def LOGP(FILE_IN): #FLOAT - GETS LOGP OF SINGLE SMILES CODE IN GIVEN FILE
    LOG_P=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    LOG_PA=LOG_P.communicate()[0].splitlines()
    for i in LOG_PA:
        if i[0:4]=="logP":
            LOG_PB=float(i.split()[1])
    return LOG_PB

def MW(FILE_IN): #FLOAT - GETS MW OF SINGLE SMILES CODE IN GIVEN FILE
    MWA=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    MWB=MWA.communicate()[0].splitlines()
    for i in MWB:
        if i[0:10]=="mol_weight":
            MWB=float(i.split()[1])
    return MWB


#PDB

def UNIPROT_SEQ(UNIPROT): #STRING - Feed UNIPROT ID and returns FASTA sequence
    import urllib2, urllib, time
    SEQ=""
    while True: #To continue fetching until it succeeds 
        try:
            url=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta"%UNIPROT).readlines()
            for i in url[1:]:
                SEQ=SEQ+i[:-1]
            break
        except urllib2.HTTPError, detail: #If HTTPError 503 is caught then it will wait 2 seconds before trying again 
            if detail.errno==503:
                time.sleep(2)
                continue
            else:
                raise
    return SEQ

def UNIPROT_CHAIN_LIMITS(UNIPROT_ID): #Given a uniprot, it returns the limits of the mature protein numbering
    import urllib2
    from Bio import SwissProt
    
    PAGE=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.txt"%UNIPROT_ID)

    PARSED_PAGE=SwissProt.parse(PAGE)
    for record in PARSED_PAGE:
        CHAIN_VALUES=[]
        for feature in record.features:
            if feature[0]=="CHAIN":
                CHAIN_VALUES=CHAIN_VALUES+[str(feature[1]), str(feature[2])]

        if any(X.isdigit()==False for X in CHAIN_VALUES) or not CHAIN_VALUES:
            CHAIN_START=1
            CHAIN_END=record.sequence_length
        else:
            CHAIN_START=min(int(X) for X in CHAIN_VALUES)
            CHAIN_END=max(int(X) for X in CHAIN_VALUES)                
    
    return[CHAIN_START, CHAIN_END]  

def GO_DICT_LIST(GO_PRE_FILE, OUT_FILE_DIC): #FILE - First file should be in the format "GO, gene ontology, GO(from is_a)
                                        #Outputs FILE of lines GO(is_a) " " GO
                                        #WARNING - MAKE SURE TO MAKE UNIQ
    FILE=open(GO_PRE_FILE)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split())    
    
    #Get rid of obsolete GOs 
    FILEA=filter(lambda x: len(x)>2,FILEA)
        
    #Get is_a TAGS:
    TAGS=[]
    for j in FILEA:
        TAGS.extend(j[2:])
    
    #MAKE DICT:
    FILEB=open(OUT_FILE_DIC,"w")
    for h in TAGS:
        for g in FILEA:
            if h in g[2:]:
                FILEB.write(h+" "+g[0]+"\n")

    
    FILEB.close()

#TEST SPEED OF DICT

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

def END_GO(GO_DICT): #Takes GO dictionary from MAKE_DICT_LIST and returns a list of terms that have no children
    DICT_VALUES=[]
    for i in GO_DICT.values():
        for j in i:
            DICT_VALUES.append(j)

    END_GO=[]
    for i in DICT_VALUES:
        if GO_DICT.has_key(i)==False:
            END_GO.append(i)

    return END_GO

def CLUSTALO_ALIGNIO(FILE_IN, FILE_ALN, FILE_OUT): # DEPRECATED!! - ALN FILE and PARSED ALN FILE - Takes in sequences saved to FILE_IN
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

def EMBOSS_ALIGNIO(FASTA1, FASTA2, FILE_OUT, LOGS):    
    import subprocess
    from Bio import AlignIO
    
    #Do global alignment
    subprocess.check_output(["/usr/local/bin/needle", FASTA1, FASTA2, 
                      "gapopen=11", "gapextend=1", "%s/EA1"%LOGS])
    
    #Convert to stockholm format for easier parsing
    FILE_IN=open("%s/EA1"%LOGS)
    ALIGNIO=AlignIO.parse(FILE_IN, "emboss")
    OUTPUT=open("%s/EA2"%LOGS, "w")
    AlignIO.write(ALIGNIO, OUTPUT, "stockholm")
    
    FILE_IN.close()
    OUTPUT.close()
    
    #Get rid of "#" and "/" containing lines
    PARSED_IN=open("%s/EA2"%LOGS)
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

def JACCARD(set1, set2): #Feed 2 sets of elements and returns float coefficient 
    result=float(len(set1&set2))/float(len(set1|set2))  
    return result

PDB_FUNCTIONS=[]

def GET_PDB(PDB_NAME,FOLDER): #PDB FILE - NEEDS Bio.PDB import*, downloads PDB_NAME to specified FOLDER
    from Bio.PDB import *
    pdb1=PDBList()
    pdb1.retrieve_pdb_file(PDB_NAME,pdir=FOLDER)
    return

def PDB_SEQ(PDB_FILE, CHAIN="A"): #STRING + LIST(int) - Give PDB file and chain and returns sequence as string and residue 
                                    #numbering as a list of integers
    #AMINO ACID DICTIONARY
    AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     #modified amino acids - Non immediate amino acid source will have value "X"
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'X', 
     'CME':'C', 'CSX':'C', 'CSO':'C', 'CSS':'C', 'SEP':'S', 'TPO':'T', 'LLP':'K',
     'SUI':'X', "NEP":"H", "FTR":"W", "TRQ":"W", "TQQ":"W", "TTQ":"W", "1TQ":"W", 
     "0AF":"W", "TRW":"W", "TOQ":"W", "OCS":"C", "CSE":"C", "TPQ":"Y", "HIP":"H",
     "CSD":"C", "PYR":"X", "PCA":"E", "SEC":"A", "MTY":"Y", "OXX":"D", "PTR":"Y",
     "PHD":"D", "MHS":"H", "AGM":"R", "MGN":"Q", "GL3":"G", "SMC":"C", "PN2":"X",
     "MLY":"K", "DM0":"K", "CAS":"C", "CSB":"C", "SEB":"S", "IAS":"X", "ACE":"X",
     "CAF":"C", "SDP":"S", "DDZ":"A", "CSR":"C", "FGL":"G", "CSW":"C"}
    if CHAIN=="":
        CHAIN="A"
    
        
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

def PDB_RES_LIG(PDB, SEQ, LOGS, CHAIN="A",D=5, IONS="TRUE"): #LIST - FOR CLASS GS -Returns numbers mapping to fed sequence of residues
                                                                #that are within D amstrongs of ligands
                                                                #PDB is ID of interest, SEQ is sequence obtained from alignment
                                                                #parsed through CLUSTALO_ALIGNIO() or any other sequence, LOGS is the LOG folder, D is
                                                                #restrictive distance, IONS=TRUE would filter out all ligands of
                                                                #length < 3  #FIRST TIME USED TO GET GOLD STANDARD IN CLASS GO PROJECT
    
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

GO_FUNCTIONS=[]

def GO_ENRICHMENT(GENE_LIST_FILE, BACKGROUND_FILE, FILE_OUT,GENE_ONTOLOGY="Process" ,p_cutoff=0.01): #GENE_LIST_FILE can be a single one 
                                                                #gene-per-line file or for batch jobs can be zipped files of this format, 
                                                                #BACKGROUN_FILE is a single gene-per-line file, FILE_OUT is where results
                                                                #will be stored, GENE_ONTOLOGY can either be: "Process", "Function" or 
                                                                #"Component", p_cutoff is the p-value cutoff for significant shared GO terms
                                                                #search
                                                                #IMPORTANT - FILES MUST HAVE FULL PATH - FATAL ERRORS MAY OCCUR IF FULL
                                                                #PATH IS NOT CORRECTLY WRITTEN                                                                
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
                                                                                    
    #Go to site
    DRIVER=webdriver.Firefox()
    DRIVER.set_window_size(600,400)    
    DRIVER.get("http://go.princeton.edu/cgi-bin/GOTermFinder")
    assert "Gene Ontology" in DRIVER.title
    
    #Introduce list of genes
    UPLOAD_GENE=DRIVER.find_element_by_name("uploadFile") #could have chosen to find it by other means, chose by_name, as long as finding
                                                            #parameter was within the "input bracket" in html obtained by page_source
    UPLOAD_GENE.send_keys(GENE_LIST_FILE)

    #Choose and ontology:
    GO=DRIVER.find_elements_by_tag_name("label") #In this instance we choose "label" because it is a checkbox, similar to select
    for j in GO:
        if j.text.lower()==GENE_ONTOLOGY.lower():j.click()
    
    #Choose annotation
    ANNOTATIONS=DRIVER.find_elements_by_tag_name("option") #again could have chosen other, but since "option" is a parameter within each 
                                                            #element of drop down menu chose it instead 
    ORGANISMS=[]
    X=1
    for organism in ANNOTATIONS: 
        ORGANISMS.append(str(X)+"#"+organism.text)
        X=X+1
    for i in ORGANISMS: print i.split("#")[0], i.split("#")[1]
    
    USER_SELECTION2=raw_input("Choose annotation by entering number key and enter:")
    
    SELECTION_ORGANISM=filter(lambda x:x.split("#")[0]==USER_SELECTION2, ORGANISMS)
    SELECTION_ORGANISM=SELECTION_ORGANISM[0].split("#")[1]
    print SELECTION_ORGANISM

    for option in ANNOTATIONS:
        if option.text==SELECTION_ORGANISM: option.click()
    
    #Upload background
    BACKGROUND=DRIVER.find_element_by_name("uploadBackground")
    BACKGROUND.send_keys(BACKGROUND_FILE)
    
    #Set p-value cut-off - 0.01 is used as default
    P_VALUE=DRIVER.find_element_by_name("pValueCutoff")
    P_VALUE.send_keys(str(p_cutoff))
    
    #SUBMIT
    SUBMIT=DRIVER.find_element_by_name("searchGOTermsButton")
    SUBMIT.click()
    
    #Check for completion and mistakes in input
    Z=0
    while Z==0:
        URL=DRIVER.current_url
        CURRENT_PAGE=[]
        for s in DRIVER.page_source.splitlines():
            CURRENT_PAGE.append(s)
        
        if "JobWatch" in URL: 
            if len(CURRENT_PAGE)==16:#WARNING - This is extremely dependent on the output format
                                     #so test with visuals every once in a while
                Z=2        
            elif len(CURRENT_PAGE)==17:
                Z=0
                URL=DRIVER.current_url
        else:
            Z=1

    #If no problems with run
    print Z
    if Z==1:
        #GET OUTPUT:
        WEB_RESULT_LINK=DRIVER.find_element_by_link_text("Tab-delimited")#This time we find by link, if you look at the HTTP it is difficult to classify
                                                                        #the link by a parameter, that's how you know it could be a link
        WEB_RESULT_LINK.click()

        DRIVER.implicitly_wait(10)        
        
        #SAVE
        WEB_RESULT=[]
        for string in DRIVER.page_source.splitlines():
            WEB_RESULT.append(string)
        WEB_RESULT=WEB_RESULT[1:-1]
        
        OUTPUT=open(FILE_OUT,"w")
        for line in WEB_RESULT:
                OUTPUT.write(line+"\n")
                
        print "File saved to %s"%FILE_OUT
    else:
        print "No known genes found. Check that genes in list are annotated with GO terms and that GO association (organism) is correct"
    
    #Close browser
    DRIVER.quit() 

METABOLITE_COMPARISSON=[]

def TC_RANK(QUERY, LIST, OUTPUT, FP="FP2"): #QUERY = Single smiles code(no name),
                                            #LIST = To compare against, in the format "SMILES    NAME" per line
                                            #OUTPUT file
    import subprocess
    
    TERMINAL_RESULT=subprocess.Popen(["-c", "/usr/local/bin/babel %s %s -ofpt -xf%s"%(QUERY, LIST, FP)], stdout=subprocess.PIPE, shell=True)
    TERMINAL_RESULT=TERMINAL_RESULT.communicate()[0].splitlines()

    #Filter out troublesome lines
    TERMINAL_RESULT=filter(lambda x:len(x)>3 and x!="Possible superstructure of first mol", TERMINAL_RESULT)

    #Rank in descending order
    TC_RANK=[]
    for tc in TERMINAL_RESULT:
        TC_RANK.append([float(tc.split("=")[-1]), tc[0:(len(tc)-len(tc.split("=")[-1]))]])
    
    TC_RANK=sorted(TC_RANK, reverse=True)    
    
    #Write file
    FILE_OUT=open(OUTPUT,"w")
    for line in TC_RANK:
        FILE_OUT.write(line[1]+" "+str(line[0])+"\n")
    FILE_OUT.close()      

def TC_RANK_FOR_106 (QUERY, SUBJECT, DICT, LOGS): #Takes metabolite NAMES that are present in a personalized
                                                    #smiles dictionary, in the case of 106.py, the dictionary
                                                    #in that script will be based on the HAVE_SMILES.smi
                                                    #WARNING - KEEP IN MIND IF DICT HAS NO KEYS FOR METABOLITES
                                                    #ENTERED THEN TC=0.0
    from FUNCTIONS import TC_RANK
    
    if SUBJECT in DICT and QUERY in DICT:
        #Format files to be taken by TC_RANK
        QUERY_FILE=open(LOGS+"/TRF1.smi", "w")
        QUERY_FILE.write(DICT[QUERY])
        QUERY_FILE.close()
        
        SUBJECT_FILE=open(LOGS+"/TRF2.smi", "w")
        SUBJECT_FILE.write(DICT[SUBJECT]+"\t"+SUBJECT)
        SUBJECT_FILE.close()
        
        #Call TC_RANK
        TC_RANK(LOGS+"/TRF1.smi", LOGS+"/TRF2.smi", LOGS+"/TRF3")
        
        #Call value from file
        FILE1=open(LOGS+"/TRF3")
        FILE_IN1=FILE1.read().splitlines()
        FILE1.close()
        
        TC=float(FILE_IN1[0].split("=")[1])
    
    else:
        TC=0.0 #IF KEY IS NOT FOUND!
    
    return TC

def TC_RANK_V2(QUERY, LIST, DICT, LOGS, FILTER=1.0, FP="FP2"): #QUERY = Name of metabolite to be analyzed
                                                #LIST = List of metabolites to be analyzed against                                        
                                                #DICT = Dictionary of MET_TO_SMILES for all entries
                                                #LOGS = Logs
                                                #FILTER= higher than it will be displayed
                                                #Returns LIST of molecules names that pass filter, or empty list if none,
                                                #if no key found in dictionary provided, it will also return an 
                                                #empty list
    import subprocess
    
    #Check first if query is in dictionary
    if DICT.has_key(QUERY):
    
        #Prep files
        QUERY_FILE=open(LOGS+"/TCRV21.smi", "w")
        QUERY_FILE.write(DICT[QUERY])
        QUERY_FILE.close()
        
        LIST_FILE=open(LOGS+"/TCRV22.smi", "w")
        for record in LIST:
            if DICT.has_key(record): #Checking for keys (metabolites or key_sets) that do not actually have smiles code
                LIST_FILE.write(DICT[record]+"\t"+record+"\n")
        LIST_FILE.close()
        
        #CHECK THAT LIST_FILE  HAS AT LEAST ONE RECORD IN IT TO RUN WITH, OTHERWISE RETURN EMPTY LIST
        CHECK_LIST=subprocess.check_output(["wc", "-l", "%s/TCRV22.smi"%LOGS])
        COUNT=int(CHECK_LIST.split()[0].strip())
        
        if COUNT>0:
        
            #Run Babel
            TERMINAL_RESULT=subprocess.Popen(["-c", "/usr/local/bin/babel %s/TCRV21.smi %s/TCRV22.smi -ofpt -xf%s"%(LOGS, LOGS, FP)], 
                                             stdout=subprocess.PIPE, shell=True)
            TERMINAL_RESULT=TERMINAL_RESULT.communicate()[0].splitlines()
        
            #Filter out troublesome lines
            TERMINAL_RESULT=filter(lambda x:len(x)>3 and x!="Possible superstructure of first mol", TERMINAL_RESULT)
            
            #Rank in descending order
            RANK=[[float(X.split("=")[1].strip()), X.split("Tanimoto")[0].strip().strip(">")] for X in TERMINAL_RESULT]
            
            RANK=sorted(RANK, reverse=True)
            
            #Return molecule names that pass filter
            FILTER=[Y[1] for Y in RANK if Y[0]>=1.0]
        
        else: #IF LIST_FILE HAS NO COMPATIBLE KEYS AT ALL IT WILL RETURN AN EMPTY LIST!!!
            FILTER=[]
    
    else: #IF QUERY IS NOT IN DICTIONARY PROVIDED IT WILL RETURN AN EMPTY LIST!!!
        FILTER=[]
        
    return FILTER

MAPPING=[]

class PDB_TO_UNIPROT_CLASS: #Given a pdb list it produces uniprot results
    import urllib2
    import pickle
    
    #Load PDB_TO_UNIPROT personal dictionarys
    DICT_FILE1=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_TO_UNIPROT_ALL.pi?attredirects=0&d=1")
    DICT=pickle.load(DICT_FILE1)
    
    #THIS IS USED FOR FUTURE UPDATES
    DICT_FILE2=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_TO_UNIPROT_ALL_SUP1.pi?attredirects=0&d=1")
    DICT.update(pickle.load(DICT_FILE2))
    ################################
    
    def __init__(self,PDB_LIST):
        self.DICT=PDB_TO_UNIPROT_CLASS.DICT
        self.PDB_LIST=PDB_LIST
        
    def has_uniprot_pdbs (self): #Displays pdbs that have a uniprot in record
        self.with_uniprot=filter(lambda x: self.DICT.has_key(x.upper())==True, self.PDB_LIST)
        return self.with_uniprot    
        
    def no_uniprot_pdbs (self): #Displays pdbs that do not have a uniprot in record
        self.without_uniprot=filter(lambda x: self.DICT.has_key(x.upper())==False, self.PDB_LIST)
        return self.without_uniprot
    
    def uniprots (self): #Displays list of matching uniprots of pdb that have uniprots #KEEP IN MIND COMPLEXES WILL HAVE FORMAT
                            #*|*
        self.uniprot_list=["|".join(self.DICT[X.upper()]) for X in self.with_uniprot]
        return self.uniprot_list
    
    def uniprot_pdb(self): #Returns list of pair [pdb,uniprot] of those pds that have a uniprot on record -FORMAT AS ABOVE
        self.results={}
        for pdb in self.with_uniprot:
            self.results[pdb.upper()]="|".join(self.DICT[pdb.upper()])
        return self.results

def PDB_TO_UNIPROT(PDB_LIST): #DICT of format {pdb:[unipro1, uniprot2,...],...}
                                #Uses http://www.uniprot.org/mapping/ server to retrieve query results as dictionary.
                                #WARNING, Only use when you are sure of the validity of your pdb list, a better choice that is
                                #curated is the class PDB_TO_UNIPROT_CLASS, use this as supplementary to the class. This function
                                #may not tell you if pdbs have been replaced by other while querying by server, use with caution 
    import urllib, urllib2

    url="http://www.uniprot.org/mapping/"

    #Set up query - It has to be as a single string of elements
    QUERY=" ".join(PDB_LIST)

    #Set up values, attributes and filters 
    params={'from':"PDB_ID", 'to':"ACC", 'format':"tab", 'query':QUERY}

    #Send job
    DATA=urllib.urlencode(params) #Encode parameters and values for form
    REQUEST=urllib2.Request(url, DATA) #Form
    RESPONSE=urllib2.urlopen(REQUEST) #Get results
    PAGE=RESPONSE.read(200000) #Process Results object into readable form

    #Process Page response
    PDB_UNIPROT_LIST=[X.split("\t") for X in PAGE.splitlines()[1:]]
    
    #Build dictionary
    PDB_UNIPROT_DICT=dict((x.upper(), []) for x in zip(*PDB_UNIPROT_LIST)[0])
    for record in PDB_UNIPROT_LIST:
        PDB_UNIPROT_DICT[record[0].upper()]=PDB_UNIPROT_DICT[record[0].upper()]+[record[1]]
    
    return PDB_UNIPROT_DICT

class UNIPROT_TO_GENES_CLASS: #Takes in list UNIPROT IDs and returns hgnc gene symbols  #NEED TO FIX FOR FALSE SYNONYMS!!!!
    "loading class of uniprot-to-gene dictionary"
    
    import urllib2
    import pickle
    FILE_IN=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_UNIPROT_TO_GENES_REVIEWED.pi?attredirects=0&d=1")
    DICT=pickle.load(FILE_IN) #Loads the uniprot_to_genes dictionary
    
    def __init__(self, UNIPROT):
        self.uniprot=UNIPROT
        self.dict=UNIPROT_TO_GENES_CLASS.DICT
    
    def uniprot_gene (self): #Returns dictionary of the form {UNIPROT:[GENE1, GENE2....]} for uniprots that have a matching gene
        self.UNIPROT_DICT={}
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                self.UNIPROT_DICT[record]=self.dict[record]
        return self.UNIPROT_DICT
    
    def has_genes_uniprots(self): #Returns lists of uniprots that have genes
        TRUE_UNIPROTS=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                TRUE_UNIPROTS.append(record)
        return TRUE_UNIPROTS
    
    def get_genes(self): #Returns list of all the genes that are present in given uniprots that have keys and are in database
        ALL_AVAILABLE_GENES=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                ALL_AVAILABLE_GENES.extend(self.dict[record])
        return ALL_AVAILABLE_GENES
        
    def no_genes(self): #Returns lists of uniprots for which genes cannot be retrieved, CHECK WITH BELOW FUNCTIONS FOR WHY
        NO_UNIPROT_GENES=[]
        for record in self.uniprot:
            if self.dict.get(record)==None or len(self.dict[record])==0:
                NO_UNIPROT_GENES.append(record)
        return NO_UNIPROT_GENES
    
    def not_on_file(self): #Returns list of uniprots that are not present in list - IMPORTANT TO ADD THEM TO FILE
        ABSENT_UNIPROTS=[]
        for record in self.uniprot:
            if self.dict.get(record)==None:
                ABSENT_UNIPROTS.append(record)
        return ABSENT_UNIPROTS
    
    def empty(self): #Returns list of uniprots that are on file but have no annotated genes - #DEPRECATED
        self.EMPTY=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])==0:
                self.EMPTY.append(record)
        return self.EMPTY

class GENES_TO_UNIPROT_CLASS: #Takes list of genes and returns uniprots. Details of output in each definition
    "Implementation for the biomaRt (R) gene id to uniprot using the hsapiens_gene_ensembl"  
    
    #Load rpy2 needed libraries
    from rpy2 import robjects
    from rpy2.robjects import r
    
    #Declare r objects
    c=robjects.r["c"]
    #Load biomaRt
    r.library("biomaRt")
    #Choose biomaRt hsapiens
    ensembl=robjects.r["useMart"]("ensembl", dataset="hsapiens_gene_ensembl")
    
    def __init__(self, GENES):
        self.robjects=GENES_TO_UNIPROT_CLASS.robjects
        self.c=GENES_TO_UNIPROT_CLASS.c
        self.GENES=GENES
        self.ensembl=GENES_TO_UNIPROT_CLASS.ensembl
        #Vectorize gene list
        self.gene_vector=self.robjects.StrVector(self.GENES)
        #Call biomaRt
        RESULTS=self.robjects.r["getBM"](attributes=self.c("hgnc_symbol", "uniprot_swissprot_accession"), filters="hgnc_symbol",
                                         values=self.gene_vector, mart=self.ensembl)
        
        #Make RESULTS pythonic
        hgnc=list(RESULTS.rx2(1))
        uniprot=list(RESULTS.rx2(2))
        self.GENE_UNIPROT_PAIRS=[list(x) for x in zip(hgnc,uniprot)]
        #To filter out record that do not have a uniprot in file
        self.PRE_GENE_UNIPROT_PAIRS_FILTERED=filter(lambda x: len(x[1])>0, self.GENE_UNIPROT_PAIRS)
        
        #Further search against personal database
        LEFTOVER_GENE_SET=filter(lambda x: x not in zip(*self.PRE_GENE_UNIPROT_PAIRS_FILTERED)[0], self.GENES)
        import urllib2
        import pickle
        REVIEWED_FILE=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_UNIPROT_TO_GENES_REVIEWED.pi?attredirects=0&d=1")
        REVIEWED_DICT=pickle.load(REVIEWED_FILE)
        
        self.PLUS_GENE_UNIPROT_PAIRS=[]
        for gene in LEFTOVER_GENE_SET:
            for key,value in REVIEWED_DICT.iteritems():
                if gene in value:
                    self.PLUS_GENE_UNIPROT_PAIRS.append([gene, key])
                elif gene.lower() in value: #to account for lower case records, if already lower case and present earlier, it
                                            #would not get to this instance
                    self.PLUS_GENE_UNIPROT_PAIRS.append([gene, key])
                    
        #Add values found in personal database to biomart found values
        self.GENE_UNIPROT_PAIRS_FILTERED=self.PRE_GENE_UNIPROT_PAIRS_FILTERED+self.PLUS_GENE_UNIPROT_PAIRS
        
    def has_uniprot_genes (self): #Returns genes in query that have a uniprot match
        GENE_PAIRS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        self.found=[]
        for gene in self.GENES:
            if gene in GENE_PAIRS:
                self.found.append(gene)
        return self.found
    
    def get_uniprots (self): #Return list of all uniprots that correspond to uniprots in has_uniprot_genes
        self.UNIPROTS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[1]
        return self.UNIPROTS
    
    def no_uniprots(self): #Return genes for which no uniprot was found #CHECK BELOW FOR REASONS
        GENE_PAIRS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        self.not_found=[]
        for gene in self.GENES:
            if gene not in GENE_PAIRS:
                self.not_found.append(gene)
        return self.not_found
    
    def not_on_database(self): #Return genes for which a record is not present in the ensembl database (hgnc) or personal
                                #database, double check them
        self.not_present=[]
        for gene in self.GENES:
            if gene not in zip(*self.GENE_UNIPROT_PAIRS)[0] and gene not in zip(*self.PLUS_GENE_UNIPROT_PAIRS)[0]:
                self.not_present.append(gene)
        return self.not_present
        
    def empty_gene(self): #Return genes that are in database but have no uniprot record - BE WARY OF RESULTS, relies on
                            #output format biomaRt, depending on size of input it may give "NA" or leave result empty giving the
                            #false impression that it is in the biomaRt database when in fact it does not exist
        self.empty=[]
        GENE_PAIRS_UNFILTERED=zip(*self.GENE_UNIPROT_PAIRS)[0]
        GENE_PAIRS_FILTERED=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        for gene in self.GENES:
            if gene in GENE_PAIRS_UNFILTERED and gene not in GENE_PAIRS_FILTERED:
                self.empty.append(gene)
        return self.empty    

class BLAST_PDB_AGAINST_CLASS: #Takes PDB|X and blasts it against subject or SWISSPROT database
    
    def __init__ (self, PDBX, IDENTITY_CUTOFF, COVERAGE_CUTOFF, LOG_FOLDER):
        #Setting variables
        self.PDBX=PDBX
        self.PDB=self.PDBX.split("|")[0].upper()
        self.CHAIN=self.PDBX.split("|")[1].upper()
        self.IC=IDENTITY_CUTOFF
        self.CC=COVERAGE_CUTOFF
        self.LOG_FOLDER=LOG_FOLDER
    
    def TO_SUBJECT(self, UNIPROT_LIST): #UNIPROT_ID STRING #If uniprot list was not provided this will produce an error
        from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, UNIPROT_CHAIN_LIMITS
        import subprocess
        
        #Download PDB
        GET_PDB(self.PDB, self.LOG_FOLDER)
        
        #Get PDB_SEQUENCE and save to fasta file if it exists
        PDB_AA, PDB_NUM=PDB_SEQ(self.LOG_FOLDER+"/"+"pdb"+self.PDB.lower()+".ent", self.CHAIN)
        
        if len(PDB_AA)>0: #To account for improperly annotated or wrong chain since it would produce an empty amino acid sequence
            PDB_AA_FASTA=open(self.LOG_FOLDER+"/PDB.fasta","w")
            PDB_AA_FASTA.write(">"+self.PDBX+"\n"+PDB_AA)
            PDB_AA_FASTA.close()
            
            #Get UNIPROT_SEQ for UNIPROT in list and save all to fasta file
            UNIPROT_AA_FASTA=open(self.LOG_FOLDER+"/UNIPROT.fasta", "w")
            for UNIPROT_ID in UNIPROT_LIST:
                UNIPROT_AA=UNIPROT_SEQ(UNIPROT_ID)
                #CHAIN_LIMITS=UNIPROT_CHAIN_LIMITS(UNIPROT_ID) #TO CORRECT FOR MATURE PROTEINS #USE WHEN LOOKING FOR COVERAGE AGAINST
                                                                                                #UNIPROT (1)
                #UNIPROT_AA=UNIPROT_AA[CHAIN_LIMITS[0]: CHAIN_LIMITS[1]] (1)
                UNIPROT_AA_FASTA.write(">"+UNIPROT_ID+"\n"+UNIPROT_AA+"\n")
            UNIPROT_AA_FASTA.close()    
                
            #Do BLAST and only show highest 
            BLAST=subprocess.check_output(["-c", "blastp -query=%s/PDB.fasta -subject=%s/UNIPROT.fasta -num_alignments=1 " 
                                           "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%(self.LOG_FOLDER, self.LOG_FOLDER)], 
                                          shell=True).splitlines()
            
            #Filter by identity>60% and coverage>80 - HERE COVERAGE=ALIGNMENT LENGTH/UNIPROT LENGTH
            if len(BLAST)>0:
                RESULT=BLAST[0].split()
                IDENTITY=float(RESULT[1])
                #COVERAGE=(float(RESULT[3])/float(RESULT[4]))*100 #USE WHEN LOOKING FOR COVERAGE AGAINST UNIPROT (1)
                COVERAGE=float(RESULT[2])
                
                if IDENTITY>float(self.IC) and COVERAGE>float(self.CC):
                    self.UNIPROT=RESULT[0]
                else:
                    self.UNIPROT="None"
            else:
                self.UNIPROT="None"
        else:
            self.UNIPROT="None"
            
        #Return answer
        return self.UNIPROT
    
    def TO_ALL (self): #UNIPROT_ID STRING 
                        #Since uniprot is not expected it aligns to all swissprot database -  THIS FUNCITON DOES TWO ROUNDS OF BLAST
        from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, UNIPROT_CHAIN_LIMITS
        import subprocess
    
        #Download PDB
        GET_PDB(self.PDB, self.LOG_FOLDER)
        
        #Get PDB_SEQUENCE and save to fasta file if it exists
        PDB_AA, PDB_NUM=PDB_SEQ(self.LOG_FOLDER+"/"+"pdb"+self.PDB.lower()+".ent", self.CHAIN)
        
        if len(PDB_AA)>0: #To account for improperly annotated or wrong chain since it would produce an empty amino acid sequence        
            PDB_AA_FASTA=open(self.LOG_FOLDER+"/PDB_ALL.fasta","w")
            PDB_AA_FASTA.write(">"+self.PDBX+"\n"+PDB_AA)
            PDB_AA_FASTA.close()
            
            #IMPORTANT - Do BLAST to get HIGHEST CANDIDATE first 
            BLAST2=subprocess.check_output(["-c", "blastp -query=%s/PDB_ALL.fasta -db=swissprot -num_alignments=1 " 
                                           "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%self.LOG_FOLDER], 
                                          shell=True).splitlines()
            
            ##########USE WHEN DOING COVERAGE AGAIN UNIPROT (1)###########                                   
            #Second blast on candidates if any
            """
            if len(BLAST)>0:
                CANDIDATE_UNIPROT=BLAST[0].split()[0].split("|")[3].split(".")[0]
                
                #Get mature sequence of candidate uniprot
                UNIPROT_CANDIDATE_AA=UNIPROT_SEQ(CANDIDATE_UNIPROT)
                CANDIDATE_CHAIN_LIMITS=UNIPROT_CHAIN_LIMITS(CANDIDATE_UNIPROT)
                UNIPROT_CANDIDATE_AA=UNIPROT_CANDIDATE_AA[CANDIDATE_CHAIN_LIMITS[0]: CANDIDATE_CHAIN_LIMITS[1]]
                
                #Save to fasta:
                UNIPROT_CANDIDATE_AA_FASTA=open(self.LOG_FOLDER+"/UNIPROT_CANDIDATE.fasta", "w")
                UNIPROT_CANDIDATE_AA_FASTA.write(">"+CANDIDATE_UNIPROT+"\n"+UNIPROT_CANDIDATE_AA+"\n")
                UNIPROT_CANDIDATE_AA_FASTA.close()
                
                #Blast
                BLAST2=subprocess.check_output(["-c", "blastp -query=%s/PDB_ALL.fasta -subject=%s/UNIPROT_CANDIDATE.fasta -num_alignments=1 "
                                                 "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%(self.LOG_FOLDER, self.LOG_FOLDER)],
                                                shell=True).splitlines() 
            
            """
            ################################################################
            
            #Filter by identity>60% and coverage>80 - HERE COVERAGE=ALIGNMENT LENGTH/UNIPROT LENGTH
            if len(BLAST2)>0:
                RESULT=BLAST2[0].split()
                IDENTITY=float(RESULT[1])
                #COVERAGE=(float(RESULT[3])/float(RESULT[4]))*100 #USE WHEN DOING UNIPROT COVERAGE (1)
                COVERAGE=float(RESULT[2])
                print COVERAGE
                
                if IDENTITY>float(self.IC) and COVERAGE>float(self.CC):
                    #self.UNIPROT=RESULT[0] #USE WHEN DOING COVERAGE AGAINST UNIPROT (1)
                    self.UNIPROT=RESULT[0].split("|")[3].split(".")[0]
                else:
                    self.UNIPROT="None"
            else:
                self.UNIPROT="None"
            
            """#USE WHEN DOING COVERAGE AGAINST UNIPROT (1)
            else:
                self.UNIPROT="None"
            """
        else:
            self.UNIPROT="None"
        #Return answer
        return self.UNIPROT 

class PDB_CHAIN_TO_UNIPROT_CLASS: #Given a list of PDB|X returns Uniprots corresponding to it
    
    def __init__ (self, PDB_LIST):
        #Load modules
        import subprocess
        from FUNCTIONS import PDB_TO_UNIPROT_CLASS, BLAST_PDB_AGAINST_CLASS
        
        #Load list
        self.PDB_LIST=[X.upper() for X in PDB_LIST]
        
        #Check for pdbs using class
        CLASS1=PDB_TO_UNIPROT_CLASS([X.split("|")[0] for X in self.PDB_LIST])
        CLASS_HAS=CLASS1.has_uniprot_pdbs()
        CLASS_HASNT=CLASS1.no_uniprot_pdbs()

        #For those pdbs that have uniprots, targeted blast function 
        self.PDB_LIST_HAS=[X for X in self.PDB_LIST if X.split("|")[0] in CLASS_HAS]
        self.PDB_CHAIN_UNIPROT={} #Dictionary to store those PDB|X that do have a uniprot
        self.PDB_CHAIN_NO_UNIPROT={} #Dictionary to store those PDB|X that do not have a uniprot, value="None"
        
        for record in self.PDB_LIST_HAS:
            CLASS2=BLAST_PDB_AGAINST_CLASS(record, 60, 80, "LOGS")
            UNIPROT_CANDIDATES=CLASS1.uniprot_pdb()[record.split("|")[0]].split("|") #Obtained from PDB_TO_UNIPROT Module
            CLASS2_RESULT=CLASS2.TO_SUBJECT(UNIPROT_CANDIDATES)

            if CLASS2_RESULT!="None":
                self.PDB_CHAIN_UNIPROT[record]=CLASS2_RESULT
            else:
                self.PDB_CHAIN_NO_UNIPROT[record]="None"
                
        #For those pdbs that do not have uniprot total blast function
        self.PDB_LIST_HASNT=[X for X in self.PDB_LIST if X.split("|")[0] in CLASS_HASNT]
        
        for record in self.PDB_LIST_HASNT:
            CLASS3=BLAST_PDB_AGAINST_CLASS(record, 60, 80, "LOGS")
            CLASS3_RESULT=CLASS3.TO_ALL()
            
            if CLASS3_RESULT!="None":
                self.PDB_CHAIN_UNIPROT[record]=CLASS3_RESULT
            else:
                self.PDB_CHAIN_NO_UNIPROT[record]="None"
    
    def with_uniprots (self): #Returns dictionary of PDB|X:uniprot that do have single qualifier uniprot
        return self.PDB_CHAIN_UNIPROT
        
    def without_uniprots (self): #Returns dictionary of PDB|X:"None" that do not have a uniprot
        return self.PDB_CHAIN_NO_UNIPROT

def PDB_RES_TO_UNIPROT_RES (PDB_DICT, LOG_FOLDER): #THROUGH *REPLACEMENT METHOD
    from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, EMBOSS_ALIGNIO, PDB_CHAIN_TO_UNIPROT_CLASS
    import pickle
    import urllib2
    import itertools
    
    #Import Biolip-based object that maps pdb-chain to uniprot #FORMAT IS {pdb|chain:uniprot, ......}
    #DICT_FILE1=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_CHAIN_TO_UNIPROT.pi?attredirects=0&d=1")
    #DICT_PDB_CHAIN_TO_UNIPROT=pickle.load(DICT_FILE1)

    UNIPROT_DICT={}
    print PDB_DICT
    for record in PDB_DICT.iterkeys():
        print record
        #Parse each record
        PDB=record.split("|")[0].upper()
        CHAIN=record.split("|")[1].upper()
        RES_LIST=PDB_DICT[record]
        
        #Get pdb file - Will download as pdb****.ent where **** is the 4 digit pdb identifier
        PDB_FILE=GET_PDB(PDB.upper(), LOG_FOLDER)
        
        #Get sequence and numbering from file
        PDB_AA, PDB_NUM=PDB_SEQ(LOG_FOLDER+"/pdb%s.ent"%PDB.lower(), CHAIN.upper())
        
        #Get Uniprot using PDB_CHAIN_TO_UNIPROT_CLASS
        UNIPROT=PDB_CHAIN_TO_UNIPROT_CLASS([PDB.upper()+"|"+CHAIN.upper()]) #Takes LIST of PDB|X, in this case just one in list
        WITH_UNIPROT=UNIPROT.with_uniprots() 
        WITHOUT_UNIPROT=UNIPROT.without_uniprots()#Those without found uniprot will return {PDB|X:"None"}    
        
        #Check if uniprot is found and continue
        if len(WITH_UNIPROT)>0:
            
            #Get Uniprot sequence
            UNIPROT_AA=UNIPROT_SEQ(WITH_UNIPROT[record.upper()])
            
            #Do asterisk replacement method - THIS WILL DO A EMBOSS PER RESIDUE, IF EVER THERE IS A BETTER WAY TO DO THIS, IMPLEMENT HERE
            UNIPROT_RES=[]
            for res in RES_LIST:
                if res in PDB_NUM: #IMPORTANT - To account for the fact that there are false pdb numbering in CSA tables                
                    #Replace with asterisk at PDB residue number
                    PDB_AST_IND=PDB_NUM.index(res)
                    PDB_AST_AA=PDB_AA[:PDB_AST_IND]+"*"+PDB_AA[PDB_AST_IND+1:]
                    
                    #Align pdb and uniprot sequences
                    PDB_FASTA=open(LOG_FOLDER+"/PRTUR1.fasta", "w") #Write sequences PER file for EMBOSS processing
                    UNIPROT_FASTA=open(LOG_FOLDER+"/PRTUR2.fasta", "w")
                    PDB_FASTA.write(">"+PDB+"\n"+PDB_AST_AA)
                    UNIPROT_FASTA.write(">"+WITH_UNIPROT[record.upper()]+"\n"+UNIPROT_AA)
                    PDB_FASTA.close()
                    UNIPROT_FASTA.close()    
                    EMBOSS_ALIGNIO(LOG_FOLDER+"/PRTUR1.fasta", LOG_FOLDER+"/PRTUR2.fasta", LOG_FOLDER+"/LOG3", LOG_FOLDER)
                    
                    #Get uniprot residue number from file
                    FILE_IN1=open(LOG_FOLDER+"/LOG3")
                    ALN_SEQUENCES=FILE_IN1.read().splitlines()
                    PDB_AL=ALN_SEQUENCES[0].split()[1]
                    UNI_AL=ALN_SEQUENCES[1].split()[1]
                    FILE_IN1.close()
                    
                    AST_COUNT=len([x for x in itertools.takewhile(lambda x: x!="*", PDB_AL)])
                    UNI_FIX=UNI_AL[:AST_COUNT].count("-") #In case of missmaches in uniprot sequence
                    UNI_COUNT=AST_COUNT+1-UNI_FIX
                    
                    UNIPROT_RES.append([UNI_COUNT,UNI_AL[AST_COUNT]])
        
            UNIPROT_DICT[record.upper()+"="+WITH_UNIPROT[record.upper()]]=UNIPROT_RES
        
        #If no uniprot found for PDB|X, will update as {PDB|X:"None"}
        else:
            UNIPROT_DICT.update(WITHOUT_UNIPROT)
        
    return UNIPROT_DICT 
   
class CSA_TO_UNIPROT_CSA_CLASS: #Dictionary Object - Make sure server is running an you have imported SQL dump file into
                                        #a database in MYSQL 
    #Load functions
    import MySQLdb as MQ
    
    def __init__ (self, SQL_DATABASE, HOST, USER, PASSWORD, LOG_FOLDER):
        from FUNCTIONS import unique, PDB_RES_TO_UNIPROT_RES        
        self.unique=unique
        self.PDB_RES_TO_UNIPROT_RES=PDB_RES_TO_UNIPROT_RES
        self.SQL_DATABASE=SQL_DATABASE
        self.HOST=HOST
        self.USER=USER
        self.PASSWORD=PASSWORD
        
        #Create connection object that represents database
        conn_obj=CSA_TO_UNIPROT_CSA_CLASS.MQ.connect(host=self.HOST,user=self.USER, passwd=self.PASSWORD, db=self.SQL_DATABASE)
        cursor=conn_obj.cursor()
        
        #Query Catalytic Residues table for PDB_ID, chain, PDB_residue_number, uniprot_residue_number, residue type 
        cursor.execute("SELECT PDBID, CHAIN, NUMBER TYPE FROM CSA_CATALYTICRESIDUES")
        PRE_RESULTS=[list(x) for x in cursor.fetchall()]
        PRE_RESULTS=self.unique(PRE_RESULTS)
        print PRE_RESULTS
        
        #Format for PDB_RES_TO_UNIPROT_RES function
        PDB_DICT=dict((x[0].upper()+"|"+x[1],[]) for x in PRE_RESULTS)
        print PDB_DICT
        for record in PRE_RESULTS:
            PDB_DICT[record[0].upper()+"|"+record[1].upper()]=PDB_DICT[record[0].upper()+"|"+record[1].upper()]+[int(record[2])]
        print PDB_DICT
        
        #Call function
        self.PDB_UNIPROT_RESULTS=self.PDB_RES_TO_UNIPROT_RES(PDB_DICT, LOG_FOLDER)

        #Find those that have UNIPROT-RES answer (match)
        self.WITH_ALL=dict((X, self.PDB_UNIPROT_RESULTS[X]) for X in self.PDB_UNIPROT_RESULTS.iterkeys() 
                           if self.PDB_UNIPROT_RESULTS[X]!="None" and len(self.PDB_UNIPROT_RESULTS[X])!=0) 
            
    def WITH_MATCH_ALL (self): #Returns dictionary in the form {"PDB|CHAIN=UNIPROT":[csa_residues]} if PDB had a match
        return self.WITH_ALL
        
    def WITH_MATCH_UNIPROT (self): #Returns uniprot dictionary of those PDB that had a match in the form {"UNIPROT""[csa_residues]}
        self.WITH_ALL_UNIPROT=dict((X.split("=")[1], self.WITH_ALL[X]) for X in self.WITH_ALL.iterkeys())
        
        return self.WITH_ALL_UNIPROT
        
    def WITHOUT_MATCH (self): #Returns dictionary in the form {"PDB|X":"None"} for those PDB that did not match a uniprot
        self.WITH_NONE=dict((X, self.PDB_UNIPROT_RESULTS[X]) for X in self.PDB_UNIPROT_RESULTS.iterkeys() 
                           if self.PDB_UNIPROT_RESULTS[X]=="None" or len(self.PDB_UNIPROT_RESULTS[X])==0)          
        return self.WITH_NONE