#CANCER_NETWORK_MAPPING_3 - DONE IN NETWORK WITH PRE-JULY DATA
#LOAD CANCER DATA TO NETWORK FOR ILLUSTRATION
import networkx as NX
import matplotlib.pyplot as plt
import pickle
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
        
#Get proteins only
PROTEINS=[]
for protein in METABO_EN.nodes():
    if METABO_EN.in_degree(protein)>0:
        PROTEINS.append(protein)
print PROTEINS
print "METABOLITES", len(METABOLITES)

#SINCE I HAVE SEPARATED METABOLITES AND PROTEINS, I HAVE TO CONVERT TO UNDIRECTED GRAPH TO GET NEIGHBORS NOW
METABO_EN=METABO_EN.to_undirected()

#Add cancer attributes to proteins
#Load cancer dict
FILE_IN1=open("DATABASES/OBJECTS/DICT_CANCER_TO_UNIPROT_GENE.pi") #HAS THE FORMAT {CANCER:[UNIPROT1|GENE, UNIPROT2|GENE...]}
CANCER_DICT=pickle.load(FILE_IN1)

#Make {UNIPROT:GENE} dictionary
UNIPROT_GENE_DICT={}
for cancer in CANCER_DICT.itervalues():
    for record in cancer:
        UNIPROT_GENE_DICT[record.split("|")[0]]=record.split("|")[1]

#BREAK DOWN NODES BY CANCER
OV_PROTEINS=[X for X in PROTEINS if X.split("|")[1] in [Y.split("|")[0] for Y in CANCER_DICT["OV"]]]
print "OV_PROTEINS", len(OV_PROTEINS)        

OV_METABOLITES=[]
for protein in OV_PROTEINS:
    OV_METABOLITES=OV_METABOLITES+METABO_EN.neighbors(protein)
OV_METABOLITES=unique(OV_METABOLITES)
print "OV_METABOLITES", len(OV_METABOLITES)

OV_EDGES=METABO_EN.edges(OV_PROTEINS)
OV_EDGES=unique(OV_EDGES)

OV_LABELS_PROTEINS={}
for protein in OV_PROTEINS:
    OV_LABELS_PROTEINS[protein]=UNIPROT_GENE_DICT[protein.split("|")[1]]

OV_LABELS_METABOLITES=dict((X[0],X[1]) for X in zip(OV_METABOLITES, OV_METABOLITES))

#DRAW NETWORK
plt.figure(figsize=(25,12))
pos=NX.spring_layout(METABO_EN)

NX.draw_networkx_nodes(METABO_EN, pos, nodelist=OV_METABOLITES, node_color="r", node_size=100)

NX.draw_networkx_nodes(METABO_EN, pos, nodelist=OV_PROTEINS, node_color="g", node_size=10)

NX.draw_networkx_edges(METABO_EN, pos, edgelist=OV_EDGES,  width=0.1)

NX.draw_networkx_labels(METABO_EN, pos, labels=OV_LABELS_PROTEINS, font_size=11, font_color="k")
NX.draw_networkx_labels(METABO_EN, pos, labels=OV_LABELS_METABOLITES, font_size=8, font_color="b")


plt.show()