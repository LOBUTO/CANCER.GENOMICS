# PDB Calculations 
import math

#Format i.PDB in pdb_files folder to keep "ATOM" and "HETATM" lines and eliminate waters
PDB= open ("PDB_LIST_LIGDIST5/Extracted_Lines/1b15") #W=5 for 1b15
PDBA=PDB.readlines()
        
print PDBA

PDBB=[]
for i in PDBA:    
    if i[0:4]=="ATOM" and i[21]=="A" or i[0:6]=="HETATM" and i[17:20]!="HOH":
        PDBB.append(i)

print PDBB
for i in PDBB: print i
print len(PDBA)
print len (PDBB) #list of all coordinates reduced by ~20% if water present 

W=5 #
X=5 #
Y=0 
SEP=0 #
Z=W+1+X+Y-SEP
print Z
CHAIN="A" #

#Get list for closest residue W
PDBC=[]
for i in PDBB: 
    if i[0:4]=="ATOM" and i[21]==CHAIN and int(i[22:26])==Z:
        PDBC.append(i) #list of coordinates for W residues

for i in PDBC: print i 

#Get substrate for closest residue W from PDB 
PDBD=[]
for j in PDBC: #coordinates at 6,7,8
    for i in PDBB:
        if i[0:6]=="HETATM":
            #IMPORTANT-IF GET MORE THAN ONE LIGAND PER PDB, MAY HAVE TO MODIFY M1 TO GET A SECOND "W" AND INCLUDE IN DIST
            DIST=math.sqrt((float(j[30:38])-float(i[30:38]))**2+(float(j[38:46])-float(i[38:46]))**2+(float(j[46:54])-float(i[46:54]))**2)
            if DIST<5 and i!=j: #Don't have to put "and i!=j" because I'm alredy filtering for HETATM
                PDBD.append(i[17:20])
                
print 6
for i in PDBD:print i

#Need to filter unique ligands




    


