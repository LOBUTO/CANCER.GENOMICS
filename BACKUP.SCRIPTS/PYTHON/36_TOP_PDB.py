#TOP_PDB
import subprocess
import shutil

def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

TOP_PDB=open("RAB_COL/TOP_PDB")
TOPA=[]
for i in TOP_PDB:TOPA.append(i.split())
TOPA=TOPA[1:]
print TOPA

PDB=[]
for i in TOPA:
    PDB.append(i[1][1:-1])
print PDB
print len(PDB)
PDBA=unique(PDB)
print len(PDBA)


ALL_PDB=subprocess.check_output("ls", cwd="PDB_FILES/PDB_LIST_COMPRESSED/").splitlines()

for i in PDBA:
    shutil.copy("PDB_FILES/PDB_LIST_COMPRESSED/pdb%s.ent.gz"%i[0:4],"RAB_COL/PDB")

ALL_PDBA=subprocess.check_output("ls", cwd="RAB_COL/PDB/").splitlines()
PDBB=[]
for i in ALL_PDBA:
    PDBB.append(i[3:7])
print PDBB
print len(PDBB)

for i in PDBA:
    if all(i[0:4]!=j for j in PDBB):
        print i
