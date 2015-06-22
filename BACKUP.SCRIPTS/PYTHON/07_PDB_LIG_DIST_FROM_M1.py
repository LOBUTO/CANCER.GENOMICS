#X-reference M1 to PDB_LIST by ID to get PDB_LIST_LIGDIST5, last used to get PDB_LIST_CSADIST 

import shutil #to copy files

M1=open("M1_CSADIST")
M1A=[]
for i in M1:
    M1A.append(i.split()) #makes lists out of M1, non-redundant list with redundant ID-X

for i in M1A: print i
print len(M1A)

M1B=[]
for i in M1A:
    i=i[0:1]
    if all(i!=j for j in M1B):
        M1B.append(i) # Gets unique list of ID_X, but not unique ID

for i in M1B: print i #236
print len(M1B)

M1C=[]
for i in M1B: 
    for j in i:
        j=j[0:4]
        #IMPORTANT - HAD TO ADD pdb 1HWX, since the time of capra list, PDBANK had replaced the name with pdb 3mw9, so the folder PDB_LIST has to copies of this pdb with these 2 different names
        shutil.copy("PDB_LIST/pdb%s.ent.gz" %j, "PDB_LIST_CSADIST") #Found 219 pdbs, since we had redundant IDs



    
