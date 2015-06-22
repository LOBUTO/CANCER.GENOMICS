#CLUSTALO - * METHOD TO GET LIGANDS USING SPDOs FROM M1_O_TOTAL

#Last modified for i.csadist5 sets using M1_CSADIST and /Extracted_Lines_CSD

import itertools
import subprocess
import math

#AMINO ACID DICTIONARY
AADICT= {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
     #modified amino acids
     'MSE':'M', 'SCY':'C', 'CCS':'C', 'KCX':'K', 'MDO':'ASG', 
     'CME':'C', 'CSX':'C','CSO':'C','CSS':'C', 'SEP':'S','TPO':'T','LLP':'K'}

M1_L_AST=[]

#CALL FROM M1
M1_O_NO1THG=open("M1_O_TOTAL")
M1=[]
for a in M1_O_NO1THG:
    M1.append(a.split()) #Get list from each set in M1

for i in M1:
    #SPDO IMPLEMENTATION - CALL IF spdO present, otherwise don't present ligand
    SPDO=open("capra_bioinf08_data/%s"%i[2])
    SPDOA=[]
    for l in SPDO:
        SPDOA.append(l)
    
    if len(SPDOA)>2: #Let's do conditional, if spdo residues are present, then calculate PDBLF for each spdO
        SPDOB=SPDOA[2:]
        
        PDBLG=[] #Where we store our spdO-Ligand pair
        for ii in SPDOB: #We do calculations for each spdO residue
            SPDON=int(ii)
        
            #SET UP THE PDB FILE
            if i[0][5]=="_":
                PDBNAME=i[0][0:5]+"A" #Fix for "_" in chain name in certain PDB IDs
            else:
                PDBNAME=i[0][0:6]
            
            #Get chain
            CHAIN=PDBNAME[5]
            
            #Get residues in chain of interest
            PDB=open("PDB_FILES/PDB_LIST_ALL/%s"%PDBNAME[0:4])
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
            M1SEQ=i[1]
            M1SEQ=M1SEQ[0:SPDON]+"*"+M1SEQ[SPDON+1:] #Replaced with asterisk at SPDO position
            
            #CLUSTALO ALIGNMENT 
            FILE=open("CLUSTALO/M1_X_PDB/%s" %(i[0]+i[2]), "w")
            FILE.write(i[0]+" " +M1SEQ+"\n"+PDBNAME+" "+ALNSEQ+"\n") #Store both sequences in file ready to be aligned
            FILE.close()
            
            CLUSTALO_STORED=subprocess.check_output(["CLUSTALO/clustalo",
                                                     "-i", "CLUSTALO/M1_X_PDB/%s" %(i[0]+i[2]),
                                                      "-o", "CLUSTALO/M1_X_PDB_ALN/%s" %(i[0]+i[2]), "--outfmt=clu","--force"])
            
            CLUSTALO=subprocess.check_output(["CLUSTALO/clustalo", "-i", "CLUSTALO/M1_X_PDB/%s" %(i[0]+i[2])]).splitlines()
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
                            PDBL.append("".join(q[17:20].split())) #appending only ligand name    
            
            PDBLF=[]
            for h in PDBL:
                if all(h!=m for m in PDBLF):
                    PDBLF.append(h)
            print PDBLF
        
            for jj in PDBLF: #Include every element per spdO of every pdb
                PDBLG.append(str(Z)+"-"+RESNAME[len(ASTX)]+"-" +jj)
        
        M1_L_AST.append(i[0:1]+i[2:3]+PDBLG)
        print M1_L_AST
    
    else: #If spdO residues do not exist
        PDBLF=['No spdO found']
        M1_L_AST.append(i[0:1]+i[2:3]+PDBLF)
        print M1_L_AST
        
#SAVE
print M1_L_AST
M1L=open("LIGAND_OUTPUT/M1_ALL_Z_L_SPDO","w")
for t in M1_L_AST:
    M1L.write("%s\n"%" ".join(t))

M1L.close()
M1_O_NO1THG.close()
SPDO.close()





    
    
    
    
    
