#Extracting all pdb files from subfolder in FOLDER pdb and placing them into PDB_LIST
import os
import glob
import shutil

A=os.listdir("pdb/") #function os.listdir() lists all directories inside a FOLDER
print A

B=[]
for i in A:
    i="pdb/"+i #to add the path to all 
    B.append(i)

print B

C=[]
for i in B:
    C.append(glob.glob("%s/*.*"%i)) #calls the the list of folder in path /pdb and lists the pathname of all the files in here
    
print C

#EXECUTE ONLY IF NEEDED, MAKE SURE TO CHANGE THE PATHS
for i in C:
    for j in i:
        shutil.copy(j,"PDB_LIST") #copies all the files "j" into folder "PDB_LIST"




