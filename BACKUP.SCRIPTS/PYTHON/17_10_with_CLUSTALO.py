#Z CALCULATION THROUGH CLUSTAL O - This calculation is independent of Pfam 

import subprocess
import itertools
import math

M1_X_PDAU=subprocess.check_output("ls",cwd="CLUSTALO/M1_X_PDAU").splitlines() #List files, parenthesis is enough because it is
                                                                                #just a command line function
                                                                                #splitlines() only works in subprcess so far
for i in M1_X_PDAU:print i

print 7

M1_L_O=[] #Matrix of 

for i in M1_X_PDAU: #List of file names in M1_X_PDAU
    FILE=open("CLUSTALO/M1_X_PDAU/%s"%i)
    FILEA=[]
    for j in FILE: FILEA.append(j.split()) #split each M1_X_PDAU file in to its M1.split and PDAU.split components
    print FILEA
    print 8
    
    M1_A=FILEA[0] #Get the M1_O component
    
    #Get W
    W=int(M1_A[4])
    
    #Get M1_SEP
    M1_SEQ=M1_A[3] #Get sequence corresponding to i.aln file
    M1_SEP=(int(M1_SEQ[0:W].count("-"))) #count number of "-" in the sequence of the M1_O component
    
    #CLUSTALO ALIGNMENT TO OBTAIN LEN
    P=subprocess.check_output(
                              ["CLUSTALO/clustalo", "-i", #do bracket because it is a specific program
                               "CLUSTALO/M1_X_PDAU/%s" %i] #capable of iterating through it, make sure %i is inside bracket
                              ).splitlines()

    #Get LEN_M1 and LEN_PDAU
    print P
    print len(P)
    PRE_LEN=len(P) #need to know number of elements in order to half them into M1_O and PDAU components
    
    PRE_LEN_M1=P[0:(PRE_LEN/2-1)] #Make list of M1 components
    print PRE_LEN_M1
    PRE_LEN_PDAU=P[(PRE_LEN/2):(PRE_LEN-1)] #Make list of PDAU components
    print PRE_LEN_PDAU
    PRE_LEN_M1_SEQ="".join(PRE_LEN_M1[1:]) #join the sequences in each component for calculations
    PRE_LEN_PDAU_SEQ="".join(PRE_LEN_PDAU[1:])
    
    PRE_LEN_M1_SEQ_C=[
                      k for k in itertools.takewhile #call .takewhile which iterates while arguments are true 
                      (lambda x: x=="-" or x.isdigit()==True, #will store "-" and "digits" till it runs into non-"-" or
                                                                #non-"digit"
                       PRE_LEN_M1_SEQ) #in M1_SEQ
                      ]
    LEN_M1=len(PRE_LEN_M1_SEQ_C)
    
    PRE_LEN_PDAU_SEQ_C=[
                      l for l in itertools.takewhile  
                      (lambda y: y=="-" or y.isdigit()==True, 
                       PRE_LEN_PDAU_SEQ) 
                      ]
    LEN_PDAU=len(PRE_LEN_PDAU_SEQ_C)
    
    #Get LEN
    LEN=LEN_M1-LEN_PDAU
    print W,M1_SEP,LEN
    
    #Get PRE_Z
    PRE_Z=W-M1_SEP+LEN+1
    
    #GET ID AND CHAIN TO CALL PDB FROM PDAU COMPONENT - This component would have the chain completion called in the absence
    #of "X" on "ID_X" done previously when seqres was called to match it
    ID_X=PRE_LEN_PDAU[0] #ID_X portion of list of PDAU components
    ID=ID_X[1:5]
    CHAIN=ID_X[6]
    
    #CALL PDB - Got from 10.py
    PDB=open("PDB_LIST_LIGDIST5/Extracted_Lines/%s" %ID) #open PDB used for alignment with CLUSTALO
    PDBA=PDB.readlines()
  
    PDBATOM=[] #list of atom lines and chain of interest
    for m in PDBA:
        if m[0:4]=="ATOM" and m[21]==CHAIN: #call CHAIN and ATOM elements that we need and get rid of other things
            PDBATOM.append(m)  #reduce file by 50% to get rid of non coordinate data and waters (-20% data)
         
    PDBHETATM=[] #list of HETATM lines without HOH
    for n in PDBA:
        if n[0:6]=="HETATM" and n[17:20]!="HOH":
            PDBHETATM.append(n) #Create separate list of only HETATM without water
    
    print PDBHETATM 

    #CONDITIONAL TO ACCOUNT FOR DIFFERENT START RESIDUE COUNT IN PDB FILE - Y CONDITIONAL USAGE
    #Get SEQRES1 - First residue of chain of interested listed in SEQRES
    PRE_SEQRES1=[]
    for aa in PDBA:
        if aa[0:6]=="SEQRES" and aa[11:12]==CHAIN:
            PRE_SEQRES1.append(aa) #Create list of SEQRES lines of CHAIN of interest from PDB file
    SEQRES1=PRE_SEQRES1[0][19:22] #Assuming first line is the first line of sequence in chain, 
    
    #Get PDB1 - First resname in PDB coordinates for chain of interest
    PDB1=PDBATOM[0][17:20] 
    
    #Y conditional
    if int(PDBATOM[0][22:26])!=1 and PDB1!=SEQRES1: #If one is different then the other one has to be different for correct
                                                    #numbering
        Y=0 #if seqres checks out no Y necessary
    else:
        #Get Y
        PRE_Y=PDBATOM[0]
        Y=int(PRE_Y[22:26])-1 #To account for disparity of sequence vs numbering of PDB residues
            
    #GET Z
    Z=PRE_Z+Y
    
    print PRE_Z,Y,Z
    
    #ACCOUNTING FOR SKIP NORMALIZATION - Some PDB files skip residue numbering, if this residue is before Z, we need to account
    #for it    
    ATOM0=len([i for i in itertools.takewhile(lambda x: x[0:4]!="ATOM", PDBA)]) # number of lines before "ATOM"
    NA=[]
    for chain in PDBA[ATOM0:]: #Line starting with ATOM coordinates till end
        if chain[21]==CHAIN: #Since HETATM and TER inside chain will have the chain identifier
                                #Some HETATM like MES(SELENOMETHIONINE) will be included in chains, we need to take them into 
                                #account in this list
            NA.append(chain) #Produces chain list, including ATOM and HETATM
    
    NB=[i for i in itertools.takewhile(lambda x: x[0:3]!="TER", NA)] #Takes while conditions are met, this is to get rid of
                                                                        #any HETATM that occurs after the end of the chain
                                                                        #TER

    NC=len([j for j in itertools.takewhile(lambda x: int(x[22:26])!=Z+1,NB)]) #Gets how far we should iterate through NB to find
                                                                            #residue numbering skips that would relevant to 
                                                                            #our Z counting Z+1
    
    COUNT1=int(NB[0][22:26]) #Getting the first residue number of NB, assuming NB's first line contains the first ATOM of chain
    N=0 #Our normalizing term for skips
    
    for nn in NB[1:NC]: #goes through elements till Z residue
        if int(nn[22:26])>COUNT1+1: #if there is a skip by a number greater than 1 (SKIP)
            N=N+(int(nn[22:26])-(COUNT1+1)) #N will store that difference
        COUNT1=int(nn[22:26]) #We update our residue count every time to keep up with residue numbering                                                                      
    
    #Normalized Z
    Z=Z+N    
    
    #Calculate distances using Z
    PDBC=[] #List of ATOM lines that correspond to Z
    for p in PDBATOM:
        if int(p[22:26])==int(Z):
            PDBC.append(p) #list of coordinates of Z atoms
    print PDBC
            
    PDBD=[] #List of HETATM lines that are withint 5A of Z atoms
    for g in PDBC:
        for q in PDBHETATM:
            if q[0:6]=="HETATM": #Don't actually need because they are all supposed to be HETATM, just checking
                #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
                DIST=math.sqrt((float(g[30:38])-float(q[30:38]))**2+(float(g[38:46])-float(q[38:46]))**2+(float(g[46:54])-float(q[46:54]))**2)
                if DIST<5 and q!=g: #Don't have to put "and q!=g" because I'm alredy filtering for HETATM
                    PDBD.append(q[17:20]) #appending only ligand name   
    
    #Filter for unique ligands
    PDBE=[]
    for r in PDBD:
        if all(r!=s for s in PDBE):
            PDBE.append(r)
    
    print PDBE #TO SHOW YOU ARE DOING SOMETHING 
    
    M1_L_O.append(M1_A[0:3]+ M1_A[4:5]+PDBE)##+M1_A[5:6])

#Save as M1_L_O
M1LS=open("LIGAND_OUTPUT/M1_L_O2","w")
for t in M1_L_O:
    M1LS.write("%s\n"%" ".join(t))
    
    
    