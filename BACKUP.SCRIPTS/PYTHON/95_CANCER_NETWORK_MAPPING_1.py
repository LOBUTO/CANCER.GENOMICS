#CANCER_NETWORK_MAPPING_1 - DONE IN NETWORK WITH PRE-JULY DATA

#TEST IF ALL CANCER GENES ARE PRESENT IN THE NETWORK

from FUNCTIONS import UNIPROT_TO_GENES_CLASS

#Get cancer genes
CANCER_GENES=open("DATABASES/CANCER_DATA/CANCER_GENES").read().splitlines()
print CANCER_GENES
print len(CANCER_GENES)

#Get NETWORK uniprots
NETWORK_UNIPROT=open("NETWORK/ENDOGENOUS_UNIPROT").read().splitlines()

#Call class
RESULTS=UNIPROT_TO_GENES_CLASS(NETWORK_UNIPROT)

#Check for presence of cancer genes in results
PRESENT=[]
ABSENT=[]
for gene in CANCER_GENES:
    if gene in RESULTS.get_genes():
        PRESENT.append(gene)
    else:
        ABSENT.append(gene)

print len(PRESENT), PRESENT
print len (ABSENT), ABSENT
