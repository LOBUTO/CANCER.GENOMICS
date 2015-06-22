#Extract all files from PDB_LIST_LIGDIST5 to PDB_LIST_LIGDIST5/Extracted_Lines (Adjusted from initial split extraction to //Extracted)
#Last used to get /Extracted_lines_CSD

import gzip #extract different sorts of files extensions
import glob #to get file names

################
PLL=[]
for i in glob.glob("PDB_LIST_CSADIST/*.ent.gz"):
    PLL.append(i)

for i in PLL:print i

for i in PLL:
    PLLA=gzip.open(i,"rb") #used to unzip .gz files listed in PLL (i.ent.gz) and store in PLLA
    
    #MODIFY TO OPEN AS LINES
    PLLC=PLLA.readlines()      
    
    #SAVING - MODIFY RANGE WHEN USING
    PLLB=open("PDB_LIST_CSADIST/Extracted_Lines_CSD/%s" %i[20:24], "w") #open file to be used to store each of them as ID.pdb        
    PLLB.writelines(PLLC) #write each unzipped filed to i.pdb 

PLLB.close()