#FOLLOW UP ON NETWORK ANALYSIS - THERE IS NOT ENOUGH INFO WITHOUT EDGE WEIGHTS
import networkx as NX
import matplotlib.pyplot as plt
from FUNCTIONS import unique

#GET NETWORK
METABO_EN=NX.read_adjlist("NETWORK/ENDOGENOUS.ADJ", create_using=NX.DiGraph(), encoding='utf-8')

#Get metabolites only
METABOLITES=[]
METABOLITES_ZERO=[]
for metabolite in METABO_EN.nodes():
    if METABO_EN.in_degree(metabolite)==0 and METABO_EN.out_degree(metabolite)>0:
        METABOLITES.append(metabolite)
    elif METABO_EN.in_degree(metabolite)==0 and METABO_EN.out_degree(metabolite)==0:
        METABOLITES_ZERO.append(metabolite)
print "METABOLITES:", len (METABOLITES)
print len (METABOLITES) + len (METABOLITES_ZERO)
print METABOLITES[:20]

#Get proteins only
PROTEINS=[]
for protein in METABO_EN.nodes():
    if METABO_EN.in_degree(protein)>0:
        PROTEINS.append(protein)
print "PROTEINS:",len (PROTEINS)
print "UNIQUE PROTEINS:", len (unique(PROTEINS))
print PROTEINS[:5]

#DRAW NETWORK
plt.figure(figsize=(25,12))
pos=NX.spectral_layout(METABO_EN)

NX.draw_networkx_nodes(METABO_EN, pos, nodelist=METABOLITES, node_color="r", node_size=100)
#NX.draw_networkx_nodes(METABO_EN, pos, nodelist=METABOLITES_ZERO, node_color="b", node_size=100) #EXCLUDE ZEROES ALTOGETHER
                                                                                                    #FOR NOW
NX.draw_networkx_nodes(METABO_EN, pos, nodelist=PROTEINS, node_color="g", node_size=10)
NX.draw_networkx_edges(METABO_EN, pos, width=0.2)

plt.show()
