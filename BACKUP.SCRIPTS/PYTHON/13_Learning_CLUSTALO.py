#CLUSTALO TEST

import subprocess

M1T=open("M1_NO1THG")
M1TA=[]
for i in M1T:
    M1TA.append(i.split())

SRU=open("SEQRES/pdb_seqres_AUNIQUE")
SRUA=[]
for j in SRU:
    SRUA.append(j.split())

CALN=subprocess.check_output(["CLUSTALO/clustalo", "-i", "SEQRES/pdb_seqres_AUNIQUE","-o", "SEQRES/TEST1.aln", "--outfmt=clu","--force"])
                            #clustalo subprocess opens (-i) file and saves to (-o) file. Give clustal format="clu" and can erase previous
                            #versions of file by "--force"
print CALN

