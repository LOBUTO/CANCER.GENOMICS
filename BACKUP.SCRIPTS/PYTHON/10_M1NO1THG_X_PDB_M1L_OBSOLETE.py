#Calculations to get ligand from PDB to M1 ID_X ########### May not be able to use since some Pfam from CAPRA don't match

#Get M1 Matrix
import math

M1_L=[] #M1-LIGAND MATRIX

M1=open("M1_NO1THG") #M1 matrix but without "1thg" info, there is no .ligdist5 or .csadist5 info in capra files
M1A=[]
for i in M1:
    M1A.append(i.split()) #632 lists 
for i in M1A:print i 

#Call PDBs from extracted list with M1 tags
for i in M1A:  #for every set in M1
    #Get ID
    M1B=i[0]
    ID=M1B[:-2] #ID without chain
    
    #Get CHAIN
    if M1B[5]=="_": #Some ID_X have empty chains (X="_"), so we have to equate them to "A" to call "Y" from pdb
        CHAIN="A"
    else:
        CHAIN=M1B[5] #For those that have them then they are the chain
    
    #Get X
    M1C=i[1:2]
    M1D=M1C[0] #get the string out of list element Pfam - OR - could have just called M1C[1] in the previous line
    M1E=[]
    M1F=0
    while M1D[M1F]!="-": #go through M1D until we meet "-" to get X number
        M1E.append(M1D[M1F])
        M1F=M1F+1
    X=int("".join(M1E)) #make sure you make X into an integer
    
    print i
    print CHAIN
    #Get W
    W=int(i[4]) #W is the 4th element in every list, make sure to make it an integer so that .count() can call it
    
    #Get SEP
    M1G=i[3]
    SEP=int(M1G[0:W].count("-")) #count the occurrence of "-" in the sequence till W, make it an integer
    
    #Get Y - Need to call PDB
    PDB=open("PDB_LIST_LIGDIST5/Extracted_Lines/%s"%ID) #Used ID previously extracted to call pdb file
    PDBA=PDB.readlines() #open pdb into lines for better processing
    
    PDBATOM=[]
    for k in PDBA:
        if k[0:4]=="ATOM" and k[21]==CHAIN: #call CHAIN and ATOM elements that we need and get rid of other things
            PDBATOM.append(k)  #reduce file by 50% to get rid of non coordinate data and waters (-20% data)
    
    PDBHETATM=[]
    for q in PDBA:
        if q[0:6]=="HETATM" and q[17:20]!="HOH":
            PDBHETATM.append(q) #Create separate list of only HETATM without water 
    
    #To see what we are doing
    print ID
    print PDBATOM[0]
    print X
    print W
    print SEP 
    
    
    PRE_Y=PDBATOM[0] #Under the assumption will be the first line of the chain that we want
    Y=int(filter(str.isdigit, PRE_Y[22:27]))-1 #To account for disparity in pdb file enumeration of residues
                                        #IMPORTANT, some resSeq#s in python have an icode letter after them, so we need
                                        #to extract the digit out the string, for example "1" out of "1C". This is under
                                        #the assumption that the following resSeq following the lists of resSeq with an 
                                        #icode attached to them will have the same resSeq#
                                        #IMPORTANT, function filter returns those elements in iterable for which the function is true,
                                        # in this case the function is str.isdigit() to check for TRUE in string. filter()
                                        #needs two elements so we place the string being called by str.isdigit() as its second
                                        #element
                                        #IMPORTANT, since readlines is being used instead of i.split(), then we don't have
                                        #to call 22:27 and could have called 22:26 instead, but we are calling the icode
                                        # as well to have usability of the function as a reference for future calls     
    print Y
    
    #Get Z 
    Z=W+1+X+Y-SEP #Z is the call number to extract pdb residues using i.aln information
    print Z
    
    #Calculate distances using Z
    PDBC=[]
    for h in PDBATOM: 
        if int(h[22:26])==Z:
            PDBC.append(h) #list of coordinates of ATOM for W residues with respect to Z
    
    PDBD=[]
    for g in PDBC: #coordinates at 6,7,8
        for l in PDBHETATM:
            if l[0:6]=="HETATM":
                #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
                DIST=math.sqrt((float(g[30:38])-float(l[30:38]))**2+(float(g[38:46])-float(l[38:46]))**2+(float(g[46:54])-float(l[46:54]))**2)
                if DIST<5 and l!=g: #Don't have to put "and l!=g" because I'm alredy filtering for HETATM
                    PDBD.append(l[17:20]) #appending only ligand name        
    
    #Filter for unique ligands
    PDBE=[]
    for m in PDBD:
        if all(m!=n for n in PDBE):
            PDBE.append(m)
    
    print PDBE #TO SHOW YOU ARE DOING SOMETHING 
    
    #Now add it to our matrix M1
    M1_L.append(i[0:3] + i[4:5] + PDBE)

#SAVE M1_L AS M1_L
M1LS=open("M1_L","w")
for p in M1_L:
    M1LS.write("%s\n"%" ".join(p))
            
            