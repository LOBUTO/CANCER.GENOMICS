#MAKE STRUCTURAL LISTS AND ASSIGN logP and m.w. to each ligand group - ANALYZE PATTERN LATER
import subprocess

def unique(LIST): #Uniqueness filter
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY  
    return LIST  
    
#Get List of ligand from both group ligand lists

#Make li

#Make list of ID_Xs that contain ligands
ALL=[]
ALN=open("LIGAND_OUTPUT/M1_ALL_Z_L_SPDO")
for i in ALN:
    i=i.split()
    if len(i)>2 and i[2]!="No": #Cause it is splittin "No sdpO found" and the i[2] element would be No
        ALL.append(i) #List of ligand containing groups
ALN.close()

#Make list of ligands excluding metals
LIGAND=[]
for j in ALL:
    j=j[2:]
    for h in j:
        h=h.split("-")
        if all(h[2]!=l for l in LIGAND) and len(h[2])>2:
            LIGAND.append(h[2]) #List of 140 pure unique ligands from ALL - Excluding metals
LIGAND=sorted(LIGAND)
print LIGAND
print len(LIGAND)

#Make list of unique ID_X
ID_X=[]
for i in ALL:
    ID_X.append(i[0])
ID_X=unique(ID_X) #List of unique ID_X, found 249
print ID_X
print len(ID_X)

#Make list of unique R-L per group - Includes metals
ID_X_R_L=[]
for i in ID_X:
    IDXRL=[i]
    for j in ALL:
        if i==j[0]:
            IDXRL.extend(j[2:])
    IDXRL=unique(IDXRL)        
    IDXRL=sorted(IDXRL[1:],key=lambda X: X[-3:])
    ID_X_R_L.append(IDXRL)
print ID_X_R_L

#Make list of ligand to residue per group

LIGC=[]
for i in ID_X_R_L:
    LIGA=[]
    for k in i:
        k=k.split("-")
        if len(k[2])>2:
            LIGA.append(k[2])
    LIGA=unique(LIGA) #Works
    for h in LIGA:
        LIGB=[h]
        for l in i:
            if h==l[-3:]:
                LIGB.append(l.split("-")[1]) #Works

        LIGC.append(LIGB) #To append each LIGB output before the next ligand gets called

LIGC=sorted(LIGC) 
print LIGC

LIGE=[]
for i in LIGAND:
    LIGD=[i]
    for j in LIGC:
        if i==j[0]:
            LIGD.append(j[1:])
    LIGE.append(LIGD)
print LIGE

#Make one list of log_P vs amino acid distribution and another list of mw vs amino acid distribution
SMILES=open("Components-smiles-oe.smi")
SMILESA=SMILES.readlines()
SMILES.close()

LIGLOGP=[]
LIGMW=[]
for i in LIGE:
    LIGF=i[0]
    for k in SMILESA:
        k=k.split()
        if i[0]==k[1]:
            LIGG=" ".join(k)    
    
    FILE=open("SMILES/%s.smi"%LIGF, "w")
    FILE.writelines(LIGG)
    FILE.close()
    
    OBPROP=subprocess.Popen(["/usr/local/bin/obprop","SMILES/%s.smi"%LIGF], stdout=subprocess.PIPE)
    OBPROPA=OBPROP.communicate()[0].splitlines()
    for j in OBPROPA:
        if j[0:4]=="logP":
            LIGLOGP.append(j.split("/")+i[1:])
        elif j[0:10]=="mol_weight":
            LIGMW.append(j.split("/")+i[1:])

print LIGLOGP #Works
print LIGMW #Works

LOGP= open("BABEL/logP","w")
for i in LIGLOGP:
    LOGP.write(str(i))    

MW=open("BABEL/mw", "w")
for j in LIGMW:
    MW.write(str(j))


