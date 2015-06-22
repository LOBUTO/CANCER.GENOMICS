#BUILT AND PRELIMINARY ANALYSIS OF NETWORK
def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

import networkx as NX
import matplotlib.pyplot as plt

#Create directed graph out of file
METABO_NT=NX.read_adjlist("NETWORK/ENDOGENOUS.ADJ", create_using=NX.DiGraph(), encoding='utf-8')

#Find nodes that have a degree of zero, meaning that we are looking for ligands that have no anotated binding partners
EMPTY_NODES=[]
for node in METABO_NT.nodes():
    if int(METABO_NT.in_degree(node))==0 and int(METABO_NT.out_degree(node))==0:
        EMPTY_NODES.append(node)

print len(EMPTY_NODES) #Found 8656

#Get distritbution of edge density for hubs(ligands) that have at least one edge
NODES_OUT=open("NETWORK/DEGREE_NONZEROES", "w")
HUB_LIST=[] #Store edge density per node
for node in METABO_NT.nodes():
    if METABO_NT.out_degree(node)>0:
        HUB_LIST.append(METABO_NT.out_degree(node)) #Looking for out_degree since those are the only ones that are ligand_to_protein
        NODES_OUT.write(str(METABO_NT.out_degree(node))+"\n")
        
print len(HUB_LIST)
print HUB_LIST[:50]

#Boxplot of degree distribution
"""
plt.boxplot(HUB_LIST)
plt.show()
"""

HUB_LIST=filter(lambda x: x<200, HUB_LIST)
plt.hist(HUB_LIST, log=True,bins=50)
plt.show()

#Get top degree nodes
TOP_NODES=[]
for node in METABO_NT.nodes():
    if int(METABO_NT.out_degree(node))>5 and int(METABO_NT.out_degree(node))<7:
        TOP_NODES.append(node)

print len(TOP_NODES)
for i in TOP_NODES: print i

