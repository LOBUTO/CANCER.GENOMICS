#HMDB_METABOLITES NETWORKX 1
#08/21/13
#Working with Networkx

import networkx as NX
import matplotlib.pyplot as plt

#GET NETWORK
METABO_EN=NX.read_adjlist("NETWORK/MEATBOLITE_TO_UNIPROT_HMDB_FILTERED.ADJ", create_using=NX.DiGraph(), encoding='utf-8')

#Get metabolites only
METABOLITES=[]
for metabolite in METABO_EN.nodes():
    if METABO_EN.in_degree(metabolite)==0 and METABO_EN.out_degree(metabolite)>0:
        METABOLITES.append(metabolite)
print len(METABOLITES)
print "unique", len(set(METABOLITES))

#Get proteins only
PROTEINS=[]
for protein in METABO_EN.nodes():
    if METABO_EN.in_degree(protein)>0:
        PROTEINS.append(protein)
print len(PROTEINS)
print "unique", len(set(PROTEINS))

#SINCE I HAVE SEPARATED METABOLITES AND PROTEINS, I HAVE TO CONVERT TO UNDIRECTED GRAPH TO GET NEIGHBORS NOW
METABO_EN=METABO_EN.to_undirected()

#Get components
COMPONENTS=NX.connected_component_subgraphs(METABO_EN)
X=0
for comp in COMPONENTS:
    print X, len(comp.nodes())
    X=X+1

#DRAW NETWORK
plt.figure(figsize=(25,12))
pos=NX.spring_layout(METABO_EN)

NX.draw_networkx_nodes(METABO_EN, pos, nodelist=METABOLITES, node_color="r", node_size=100)

NX.draw_networkx_nodes(METABO_EN, pos, nodelist=PROTEINS, node_color="b", node_size=10)

NX.draw_networkx_edges(METABO_EN, pos, width=0.1)

plt.show()
