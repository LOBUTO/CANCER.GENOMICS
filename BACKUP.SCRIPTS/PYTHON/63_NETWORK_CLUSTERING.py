#CLUSTER NETWORK BY DEGREE QUANTILES

import networkx as NX
import subprocess
import os

#Create directed graph
METABO_NT=NX.read_adjlist("NETWORK/ENDOGENOUS.ADJ", create_using=NX.DiGraph(), encoding="utf-8")

#Get rid of hubs that have a degree of zero
NON_ZERO_NODES=[] #As before found 2319
for node in METABO_NT.nodes():
    if int(METABO_NT.out_degree(node))>0:
        NON_ZERO_NODES.append(node)

print len(NON_ZERO_NODES)

#Break the non-zero nodes by degree clustering and save protein from hub groups into files
LOW_DEGREE, MIDDLE_DEGREE, HIGH_DEGREE, NON_HUBS=[],[],[],[]
NON_HUB_PROTEINS=open("NETWORK/NON_HUB_PROTEINS","w")
LOW_DEGREE_PROTEINS=open("NETWORK/LOW_DEGREE_PROTEINS","w")
MIDDLE_DEGREE_PROTEINS=open("NETWORK/MIDDLE_DEGREE_PROTEINS","w")
HIGH_DEGREE_PROTEINS=open("NETWORK/HIGH_DEGREE_PROTEINS","w")

for node in NON_ZERO_NODES:
    NEIGHBORS=[]
    if int(METABO_NT.out_degree(node))>2 and int(METABO_NT.out_degree(node))<9:
        LOW_DEGREE.append(node)
        NEIGHBORS.extend(METABO_NT.neighbors(node))
        for nb in NEIGHBORS:
            LOW_DEGREE_PROTEINS.write(nb+"\n")
    elif int(METABO_NT.out_degree(node))>8 and int(METABO_NT.out_degree(node))<42:
        MIDDLE_DEGREE.append(node)
        NEIGHBORS.extend(METABO_NT.neighbors(node))
        for nb in NEIGHBORS:
            MIDDLE_DEGREE_PROTEINS.write(nb+"\n")
    elif int(METABO_NT.out_degree(node))>41:
        HIGH_DEGREE.append(node)
        NEIGHBORS.extend(METABO_NT.neighbors(node))
        for nb in NEIGHBORS:
            HIGH_DEGREE_PROTEINS.write(nb+"\n")
    elif int(METABO_NT.out_degree(node))<3:
        NON_HUBS.append(node)
        NEIGHBORS.extend(METABO_NT.neighbors(node))
        for nb in NEIGHBORS:
            NON_HUB_PROTEINS.write(nb+"\n")

NON_HUB_PROTEINS.close()
LOW_DEGREE_PROTEINS.close()
MIDDLE_DEGREE_PROTEINS.close()
HIGH_DEGREE_PROTEINS.close()

print len(NON_HUBS)
print len(LOW_DEGREE)
print len(MIDDLE_DEGREE)
print len(HIGH_DEGREE)
print len(LOW_DEGREE)+len(MIDDLE_DEGREE)+len(HIGH_DEGREE)+len(NON_HUBS)

#Filter files for uniqueness and store only a list of the uniprot ids
os.chdir("NETWORK")

DEGREE_FILES=["NON_HUB_PROTEINS", "LOW_DEGREE_PROTEINS", "MIDDLE_DEGREE_PROTEINS", "HIGH_DEGREE_PROTEINS"]

for file in DEGREE_FILES:
    print file
    subprocess.call('''sort %s | uniq | awk 'BEGIN {FS="|"} ; {print $2}' > TEST'''%file, shell=True)
    subprocess.call("rm %s"%file, shell=True)
    subprocess.call("mv TEST %s"%file, shell=True)

