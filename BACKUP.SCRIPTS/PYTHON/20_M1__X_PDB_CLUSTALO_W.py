#CLUSTALO - * METHOD TO GET LIGANDS USING Ws FROM SETS THAT HAVE i.aln FILES

import itertools
import subprocess
import math

#AMINO ACID DICTIONARY
AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'ASG', 'CME':'C', 'CSX':'C','CSO':'C'} #Account for modified amino acids

M1_L_AST=[]

#CALL FROM M1
M1_O_NO1THG=open("M1_O_NO1THG")
M1=[]
for a in M1_O_NO1THG:
    M1.append(a.split()) #Get list from each set in M1

for i in M1:
    #Get W
    W=int(i[4])

    #SET UP THE PDB FILE
    if i[0][5]=="_":
        PDBNAME=i[0][0:5]+"A" #Fix for "_" in chain name in certain PDB IDs
    else:
        PDBNAME=i[0][0:6]
    
    #Get chain
    CHAIN=PDBNAME[5]
    
    #Get residues in chain of interest
    PDB=open("PDB_LIST_LIGDIST5/Extracted_Lines/%s"%PDBNAME[0:4])
    PDBA=PDB.readlines()
    
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
    
    #Make sequence out of RESNAME using aa dictionary
    ALNSEQ=""
    for e in RESNAME:
        ALNSEQ=ALNSEQ+AADICT[e] #This is supposed to be the correct PDB refseq that matches PDB numbering
    
    #DO ASTERISK REPLACEMENT IN i.aln    
    #Get i.aln sequence
    M1SEQ=i[3]
    M1SEQ=M1SEQ[0:W]+"*"+M1SEQ[W+1:] #Replaced with asterisk at W position
    
    #CLUSTALO ALIGNMENT 
    FILE=open("CLUSTALO/M1_X_*BUILT/%s" %(i[0]+i[5]), "w")
    FILE.write(i[0]+" " +M1SEQ+"\n"+PDBNAME+" "+ALNSEQ+"\n") #Store both sequences in file ready to be aligned
    FILE.close()
    
    CLUSTALO_STORED=subprocess.check_output(["CLUSTALO/clustalo",
                                             "-i", "CLUSTALO/M1_X_*BUILT/%s" %(i[0]+i[5]),
                                              "-o", "CLUSTALO/M1_X_ALN/%s" %(i[0]+i[5]), "--outfmt=clu","--force"])
    
    CLUSTALO=subprocess.check_output(["CLUSTALO/clustalo", "-i", "CLUSTALO/M1_X_*BUILT/%s" %(i[0]+i[5])]).splitlines()
    print CLUSTALO
    #Get M1 and PDB component of alingment
    CLUSTALOM1=CLUSTALO[0:(len(CLUSTALO)/2)]
    CLUSTALOPDB=CLUSTALO[(len(CLUSTALO)/2):]
    
    #CAVEAT - To account for certain annotations placed at the output of CLUSTALO
    CAVEATA=[w for w in itertools.takewhile(lambda y: y!=">%s"%i[0][0:6],CLUSTALOM1)]
    CLUSTALOM1=CLUSTALOM1[len(CAVEATA):]
    
    #Join sequences
    SEQM1=CLUSTALOM1[0:1] + "".join(CLUSTALOM1[1:]).split()
    SEQPDB=CLUSTALOPDB[0:1] +"".join(CLUSTALOPDB[1:]).split()
    print SEQM1
    #GET Z FROM ATERISK COUNT
    ASTX="".join([f for f in itertools.takewhile(lambda z: z!="*", SEQM1[1])]) #Sequence till asterisk in i.aln from M1
    print ASTX
    Z=RESNUMBER[len(ASTX)]#This works because the numbering was obtained directly from the sequence we aligned, somewhat like
                            #a dictionary so it would call the exact number of the residue it needed
        
    #CALCULATIONS AS IN 17 AND 10, NO NEED FOR Y             
    PDBHETATM=[] #list of HETATM lines without HOH
    for n in PDBA:
        if n[0:6]=="HETATM" and n[17:20]!="HOH":
            PDBHETATM.append(n) #Create separate list of only HETATM without water 
    
    PDBZ=[]
    for f in PDBD: #Use chain atoms including HETATM embedded if any
        if int(f[22:26])==int(Z):
            PDBZ.append(f) #list of coordinates of Z atoms
    
    print i
    PDBL=[] #List of HETATM lines that are withint 5A of Z atoms
    for g in PDBZ:
        for q in PDBHETATM:
            if q[0:6]=="HETATM": #Don't actually need because they are all supposed to be HETATM, just checking
                #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
                DIST=math.sqrt((float(g[30:38])-float(q[30:38]))**2+(float(g[38:46])-float(q[38:46]))**2+(float(g[46:54])-float(q[46:54]))**2)
                if DIST<5 and q!=g: #Don't have to put "and q!=g" because I'm alredy filtering for HETATM
                    PDBL.append(q[17:20]) #appending only ligand name    
    
    PDBLF=[]
    for h in PDBL:
        if all(h!=m for m in PDBLF):
            PDBLF.append(h)
    print PDBLF
    
    M1_L_AST.append(i[0:3]+i[4:5]+PDBLF+str(Z).split())

#SAVE
print M1_L_AST
M1L=open("LIGAND_OUTPUT/M1_L_ASTK","w")
for t in M1_L_AST:
    M1L.write("%s\n"%" ".join(t))

M1L.close()