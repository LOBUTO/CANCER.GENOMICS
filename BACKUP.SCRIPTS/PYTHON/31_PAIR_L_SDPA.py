#BUILDING SDPA_COMMON LIGAND DATASET TO PAIR_L_SDPA
import itertools
import subprocess
import math

def unique(LIST):
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

#AMINO ACID DICTIONARY
AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     #modified amino acids
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'ASG', 
     'CME':'C', 'CSX':'C','CSO':'C','CSS':'C', 'SEP':'S','TPO':'T','LLP':'K'}

#Import UNICODEB
UNICODEB=open("UNIPROT_PDB/UNICODEB").readlines()

PAIR=subprocess.check_output("ls", cwd="BIOPYTHON/ALIGNIO_CAPRA_PAIR/").splitlines()

ALL_PAIRS=[]
for i in PAIR:
    print i
    PAIR1=open("BIOPYTHON/ALIGNIO_CAPRA_PAIR/%s"%i).readlines()
    print PAIR1
    X=0
    #Replacement of UNIPROT ID by PDB ID in line
    while X!=len(PAIR1):
        for h in UNICODEB:
            if any(PAIR1[X].split("|")[0]==l for l in h.split()):
                print h.split()[0].lower() 
                PAIR1[X]=h.split()[0].lower()+"|"  +"".join(PAIR1[X].split("|")[2:])
        X=X+1
    print PAIR1
    
    #Call sdpL #s corresponding to pair
    SDPA=open("capra_bioinf08_data/%ssdpA"%i[:-4]).readlines()
    SDPA=filter(lambda a: a[0]!="#", SDPA)
    print SDPA
    
    #Account for the fact that not all i.sdpA files have sdpA coordinates:
    if len(SDPA)>0:
        #GET A LIGAND FOR EACH SDPA#, THEN WE COMPARE IF THEY ARE THE SAME
        
        #DO PER ID IN PAIR1 
        
        #Set up PDB name to call PDB file
        LIG_PAIR=[]
        for ii in PAIR1:
            print ii
            PDB=ii.split("|")[0]
            if len(PDB)<5:
                ID_X=PDB+"_A"
            elif PDB[5]=="_":
                ID_X=PDB[0:5]+"A"
            else:
                ID_X=PDB
            print ID_X
            
            #Get chain:
            if PDB=="1pfx" or PDB=="1aut": CHAIN="C" #MINOR MODIFICATION TO ACCOUNT FOR PDB FILE THAT DOES NOT HAVE ALL CHAINS
            elif PDB=="1c5m": CHAIN="D"
            else:CHAIN=ID_X[5]
            
            #For ID_X per spdA
            LIG_SDPA=[] # STORE LIGANDS PER ID_X
            RES_SDPA=[]
            for iii in SDPA:
                #Get integer sdpA
                SDPAN=int(iii)
            
                #Get residues in chain of interest
                PDB_FILE=open("PDB_FILES/PDB_LIST_ALL/%s"%ID_X[0:4])
                PDBA=PDB_FILE.readlines()
                
                PDBB=PDBA[len([r for r in itertools.takewhile(lambda x: x[0:4]!="ATOM", PDBA)]):]#Get coordinate atoms plus anything that follows
                
                PDBC=[]
                for b in PDBB:
                    if b[21:22]==CHAIN and b[0:6]!="ANISOU": #Some pdb files have ANISOU embedded in between ATOM counts
                        PDBC.append(b) #Get residues that are relevant to chain of interest including HETATM if present plus following HETATM
                
                PDBD=[j for j in itertools.takewhile(lambda x: x[0:3]!="TER", PDBC)] #Get residues of chain of interest including HETATM if any
                
                #Get resname and resnumber lists    
                PDBE=[]
                for c in PDBD:
                    if all((c[17:20]+c[22:26])!=d for d in PDBE):
                        PDBE.append(c[17:20]+c[22:26]) #Get unique list of resname-resnumber pair
                    
                RESNAME=[]
                RESNUMBER=[]
                for d in PDBE:
                    RESNAME.append(d[0:3]) #List of resnames
                    RESNUMBER.append(int(d[3:])) #List of resnumbers                    
                print RESNAME
                                    
                #Make sequence out of RESNAME using aa dictionary
                ALNSEQ=""
                for e in RESNAME:
                    ALNSEQ=ALNSEQ+AADICT[e] #This is supposed to be the correct PDB refseq that matches PDB numbering
                print ALNSEQ    
                
                #DO ASTERISK REPLACEMENT IN PAIR(ii)
                M1SEQ=ii.split()[1]
                M1SEQ=M1SEQ[0:SDPAN]+"*"+M1SEQ[SDPAN+1:]
                
                #CLUSTALO ALIGNMENT
                FILE=open("CLUSTALO/PAIR_X_PDB/%s"%(str(PDB)+"_"+str(SDPAN)+"_"+i), "w")
                FILE.write(PDB+" "+M1SEQ+"\n"+ID_X+" "+ALNSEQ+"\n")
                FILE.close()
                
                CLUSTALO_STORED=subprocess.check_output(["CLUSTALO/clustalo",
                                                  "-i", "CLUSTALO/PAIR_X_PDB/%s"%(PDB+"_"+str(SDPAN)+"_"+i),
                                                  "-o", "CLUSTALO/PAIR_X_PDB_ALN/%s"%(PDB+"_"+str(SDPAN)+"_"+i),
                                                   "--outfmt=clu", "--force"])
                
                CLUSTALO=subprocess.check_output(["CLUSTALO/clustalo", 
                                                  "-i", "CLUSTALO/PAIR_X_PDB/%s"%(PDB+"_"+str(SDPAN)+"_"+i)]).splitlines()
                print CLUSTALO
                
                #Get M1(asterisk composed at M1SEQ) and PDB component of alignment
                CLUSTALOM1=CLUSTALO[0:(len(CLUSTALO)/2)]
                CLUSTALOPDB=CLUSTALO[(len(CLUSTALO)/2):]                
                
                #CAVEAT - To account for certain annotations placed at the output of CLUSTALO
                CAVEATA=[w for w in itertools.takewhile(lambda y: y!=">%s"%PDB[0:6],CLUSTALOM1)]
                CLUSTALOM1=CLUSTALOM1[len(CAVEATA):]
                
                #Join sequences
                SEQM1=CLUSTALOM1[0:1] + "".join(CLUSTALOM1[1:]).split()
                SEQPDB=CLUSTALOPDB[0:1] +"".join(CLUSTALOPDB[1:]).split()
                print SEQM1                           
                
                #GET Z FROM ATERISK COUNT
                ASTX="".join([f for f in itertools.takewhile(lambda z: z!="*", SEQM1[1])]) #Sequence till asterisk in i.aln from M1
                print ASTX
                Z=RESNUMBER[len(ASTX)] - SEQPDB[:len(ASTX)].count("-") #In case sequences don't overlap perfectly              
                
                #CALCULATIONS
                PDBHETATM=[] #list of HETATM lines without HOH
                for n in PDBA:
                    if n[0:6]=="HETATM" and n[17:20]!="HOH":
                        PDBHETATM.append(n) #Create separate list of only HETATM without water    
                
                PDBZ=[]
                for f in PDBD: #Use chain atoms including HETATM embedded if any
                    if int(f[22:26])==int(Z):
                        PDBZ.append(f) #list of coordinates of Z atoms
                
                PDBL=[] #List of HETATM lines that are withint 5A of Z atoms
                for g in PDBZ:
                    for q in PDBHETATM:
                        if q[0:6]=="HETATM": #Don't actually need because they are all supposed to be HETATM, just checking
                            #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
                            DIST=math.sqrt((float(g[30:38])-float(q[30:38]))**2+(float(g[38:46])-float(q[38:46]))**2+(float(g[46:54])-float(q[46:54]))**2)
                            if DIST<5 and q!=g: #Don't have to put "and q!=g" because I'm already filtering for HETATM
                                PDBL.append("".join(q[17:20].split())) #appending only ligand name     

                PDBLF=[]
                for h in PDBL:
                    if all(h!=m for m in PDBLF) and len(h)>2: #Get rid of ions
                        PDBLF.append(h)
                print PDBLF
                
                #Get all ligands found for all sdpA#
                LIG_SDPA.extend(PDBLF)
                
                #Get all sdpA residues per ID_X
                RES_SDPA.append(str(Z)+"-"+RESNAME[len(ASTX)])
            
            #Keep common ligands found across sdpA#s
            print LIG_SDPA
            LIG_SDPA_ALL=[]
            for o in LIG_SDPA:
                if LIG_SDPA.count(o)==len(SDPA): #If ligand count in ligand list is equal to number of sdpa
                    LIG_SDPA_ALL.append(o)
            LIG_SDPA_ALL=unique(LIG_SDPA_ALL)
            
            #GET LIGAND AND RESIDUES CORRESPONDING TO IT DEPENDING ON sdpA#s
            if len(LIG_SDPA_ALL)>0:
                LIG_PAIR.append(LIG_SDPA_ALL+RES_SDPA)
            else:
                LIG_PAIR.append("No common ligand")
            print LIG_PAIR
            
        ALL_PAIRS.append(i.split()+LIG_PAIR)
    else: 
        ALL_PAIRS.append(i+" "+"No sdpA found")
        #No sdpA#s found in i.sdpA
    print ALL_PAIRS        

FILE_OUT=open("LIGAND_OUTPUT/PAIR_L_SDPA_NO2","w")
for i in ALL_PAIRS:
    FILE_OUT.write("%s\n"%str(i))  #WORKS                
                