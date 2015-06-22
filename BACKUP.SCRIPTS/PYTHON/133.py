#BACKGROND MODEL PROJECT 2
#Calculate phylogenetic distance between "Core Patients" and rest of patients
#012614

from Bio import Phylo

#Load file - Nexus rather than Newick, Newick data from R was not identifying node names correctly
BRCA_CORR_TREE=Phylo.read("DATABASES/CANCER_DATA/TCGA/FROM_R_ANALYSIS/012614_BRCA_EXP_TREE.nex", "nexus")

#Get terminal nodes - (Basically all nodes)
TERMINAL_NODES=[x.name for x in BRCA_CORR_TREE.get_terminals()]

#Load BRCA top patients numeric file
FILE_IN1=open("DATABASES/CANCER_DATA/TCGA/FROM_R_ANALYSIS/122614_BRCA_TOP_CORR_PATIENTS_NUMERIC")
BRCA_TOP_NUMERIC=[X for X in FILE_IN1.read().splitlines()]
FILE_IN1.close()
print BRCA_TOP_NUMERIC

#Get average distances from all NON_TOP nodes to TOP NODES and write to FILE
FILE_OUT1=open("DATABASES/CANCER_DATA/TCGA/012614_BRCA_NON_TOP_NUMERIC_DISTANCE_TO_TOP","w")
for node in TERMINAL_NODES:
    if node not in BRCA_TOP_NUMERIC:
        AVG_TO_TOP=sum([BRCA_CORR_TREE.distance(node,X) for X in BRCA_TOP_NUMERIC])/len(BRCA_TOP_NUMERIC)
        
        #Write to file
        FILE_OUT1.write(str(node)+"\t"+str(AVG_TO_TOP)+"\n")

FILE_OUT1.close()



