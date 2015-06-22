#CREATE LIST OF PDBs THAT MATCH UNIPROT IDs IN THE i.aln files TO PDB_FILES/PDB_UNIPROT/

import gzip
import subprocess

#Get dictionary of PDB-UNIPROT(s) codes
DICTA=[]
DICT=open("UNIPROT_PDB/UNICODEB")
for i in DICT:
    DICTA.append(i.split())
print DICTA

#Get i.aln file list
ALN=subprocess.check_output("ls", cwd="capra_bioinf08_data").splitlines()
ALNA=[]
for i in ALN:
    if i[-3:]=="aln":
        ALNA.append(i)

for i in ALNA:
    ALNB=open("capra_bioinf08_data/%s"%i)
    ALNC=ALNB.readlines()
    for j in ALNC:
        for l in DICTA:
            if any(j.split("|")[0]==k for k in l):
                PDB=gzip.open("PDB_FILES/PDB_LIST_COMPRESSED/pdb%s.ent.gz"%l[0].lower(),"rb")
                PDBA=PDB.readlines()
                PDBB=open("PDB_FILES/PDB_UNIPROT/%s"%l[0].lower(),"w")
                PDBB.writelines(PDBA)
                PDB.close()
                PDBB.close()
                
            
    