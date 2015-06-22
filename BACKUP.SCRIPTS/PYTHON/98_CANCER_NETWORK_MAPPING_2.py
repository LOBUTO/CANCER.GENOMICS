#CANCER_NETWORK_MAPPING_2 - DONE IN NETWORK WITH PRE-JULY DATA
#Check how many genes are actual present in network

import pickle
from FUNCTIONS import UNIPROT_TO_GENES_CLASS, unique

#Load CANCER_to_Gene dict
FILE_IN1=open("DATABASES/OBJECTS/DICT_CANCER_TO_GENES.pi")
CANCER_DICT=pickle.load(FILE_IN1)
FILE_IN1.close()

#Get all cancer genes
ALL_CANCER_GENES=[]
for genes in CANCER_DICT.itervalues():
    ALL_CANCER_GENES=ALL_CANCER_GENES+genes
ALL_CANCER_GENES=unique(ALL_CANCER_GENES)
print len(ALL_CANCER_GENES)

#Load NETWORK UNIPROTS
FILE_IN2=open("NETWORK/ENDOGENOUS_UNIPROT")
UNIPROTS=FILE_IN2.read().splitlines()

#Run class
UNIPROT_CLASS=UNIPROT_TO_GENES_CLASS(UNIPROTS)
NETWORK_GENES=UNIPROT_CLASS.get_genes()
NETWORK_UNIPROT_GENES=UNIPROT_CLASS.uniprot_gene()
print "uniprot in network has match", len(UNIPROT_CLASS.has_genes_uniprots())
print "total cancer genes in network", len([X for X in ALL_CANCER_GENES if X in NETWORK_GENES])

#Check per cancer and make UNIPROT_DICT
CANCER_TO_UNIPROT_DICT=dict((X, []) for X in CANCER_DICT.iterkeys())
print CANCER_TO_UNIPROT_DICT

for cancer in CANCER_DICT.iterkeys():
    TOTAL=len(CANCER_DICT[cancer])
    PRESENT=[X for X in CANCER_DICT[cancer] if X in NETWORK_GENES]
    #print cancer, len(PRESENT), TOTAL
    
    UNIPROT_GENE_LIST=[]
    for gene in PRESENT:
        for uniprot, genes in NETWORK_UNIPROT_GENES.iteritems():
            if gene in genes:
                UNIPROT_GENE_LIST.append(uniprot+"|"+gene)
    
    CANCER_TO_UNIPROT_DICT[cancer]=UNIPROT_GENE_LIST

#CHECK DICT
for cancer, record in CANCER_TO_UNIPROT_DICT.iteritems():
    print cancer, record
#for cancer, record in CANCER_TO_UNIPROT_DICT.iteritems():
#    print cancer, len(record)

#SAVE TO OBJECT
PICKLE_OUT=open("DATABASES/OBJECTS/DICT_CANCER_TO_UNIPROT_GENE.pi", "w")
pickle.dump(CANCER_TO_UNIPROT_DICT, PICKLE_OUT)
    
    