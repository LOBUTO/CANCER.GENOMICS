###ANALYZE NATURE OF METABOLITES IN NETWORK THAT ARE LEAST 90% TO HYALURONAN

from FUNCTIONS import TC_RANK, unique
import subprocess
import networkx as NX
import matplotlib.pyplot as plt
import os

def UNIPROT_GENE(UNIPROT):#STRING - Returns gene name
    import urllib,urllib2
    url=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta"%UNIPROT).readlines()
    for line in url[0].split():
        if "GN=" in line:
            GN=line.split("=")[1]
    return GN

def UNIPROT_GENE_PLUS(UNIPROT): #LIST-The difference between this and UNIPROT_GENE is that UNIPROT_GENE_PLUS returns synonim genes as well if    
                                #any and the gene name in the first entry
    import urllib, urllib2
    from Bio import SwissProt
    
    url=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.txt"%UNIPROT)
    GENES=[]
    for record in SwissProt.parse(url):
        if len(record.gene_name.split(";"))>2:
            GENES.append(record.gene_name.split(";")[0].split("=")[1])
            SYN=record.gene_name.split(";")[1].split("=")[1].split(",")
            for syno in SYN:
                GENES.append("".join(syno.split()))
        else:
            GENES.append(record.gene_name.split(";")[0].split("=")[1])
    return GENES
    
SEPARATOR=[]

#Get Network - save degree distribution of whole network to file
TO_R1=open("/Users/jzamalloa/Desktop/FOLDER/LAB/ECLIPSE_TO_R/NETWORK_DEGREE", "w")
METABO_NT=NX.read_adjlist("NETWORK/ENDOGENOUS.ADJ", create_using=NX.DiGraph(), encoding="utf-8")
for node in METABO_NT.nodes():
    if METABO_NT.in_degree(node)==0: #Cause we just want metabolites and not proteins
        print node
        TO_R1.write(str(METABO_NT.out_degree(node))+"\n")
TO_R1.close()

#Get List of Metabolites that are at least 90% similar to Hyaluronan including Hyaluronan
TC_RANK("/Users/jzamalloa/Desktop/FOLDER/LOGS/HMDB10366.smi",
        "/Users/jzamalloa/Desktop/FOLDER/LOGS/END_SMILES.smi",
        "LOGS/TC_HN")

TC_LIST=open("LOGS/TC_HN").read().splitlines()
TC_LIST=filter(lambda x: float(x.split("=")[1])>0.90 ,TC_LIST)

##List degree of these metabolites in network - Save degree distribution of HUBS_90 to file
CLEAN_TC_LIST=[]
for line in TC_LIST:
    CLEAN_TC_LIST.append("_".join(line.split("Tanimoto")[0].split())[1:])
print len(CLEAN_TC_LIST)

TO_R2=open("/Users/jzamalloa/Desktop/FOLDER/LAB/ECLIPSE_TO_R/HUBS_90_DEGREE", "w")
HUBS_90_DEGREE=[]
HUBS_90_PROTEINS=[]
for metabolite in CLEAN_TC_LIST:
    HUBS_90_DEGREE.append(METABO_NT.out_degree(metabolite))
    TO_R2.write(str(METABO_NT.out_degree(metabolite))+"\n")
    PROTEINS=METABO_NT.neighbors(metabolite)
    for prot in PROTEINS:
        HUBS_90_PROTEINS.append(prot.split("|")[1])
TO_R2.close()

print len(HUBS_90_DEGREE)
print len(HUBS_90_PROTEINS)
print HUBS_90_PROTEINS[:5]
HUBS_90_PROTEINS=unique(HUBS_90_PROTEINS)
print len(HUBS_90_PROTEINS)

#Plot degree distribution
plt.hist(HUBS_90_DEGREE)
plt.show()

#Get gene names out of these proteins, use uniprot_sprot.dat
HUBS_90_GENES=[]
for protein in HUBS_90_PROTEINS:
    PROT_GENE=UNIPROT_GENE_PLUS(protein)
    for gene in PROT_GENE:
        HUBS_90_GENES.append(protein+"$"+gene)
    
print HUBS_90_GENES[0:5]
print len(HUBS_90_GENES)

#Get list of cancer genes:
os.chdir("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/geneStats")
CANCER_GENES=subprocess.Popen(["-c", '''ls | awk 'BEGIN {FS="."} ; {print $1}' '''], 
                              stdout=subprocess.PIPE, shell=True).communicate()[0].splitlines()
    
#Compare list of "90" genes to cancer genes
CANCER_90=[]
for protein_gene in HUBS_90_GENES:
    if protein_gene.split("$")[1] in CANCER_GENES:
        CANCER_90.append(protein_gene)

print CANCER_90
print len(CANCER_90)
    
#Find out what metabolites these proteins are related to in network group of HUBS_90 METABOLITES (CLEAN_TC_LIST)
CANCER_90_METABOLITES=[]
for protein_gene in CANCER_90:
    for node in CLEAN_TC_LIST:
        if any(protein_gene.split("$")[0]==protein.split("|")[1] for protein in METABO_NT.neighbors(node)):
            CANCER_90_METABOLITES.append(node)

print len(CANCER_90_METABOLITES)
print CANCER_90_METABOLITES
CANCER_90_METABOLITES=unique(CANCER_90_METABOLITES)
print len(CANCER_90_METABOLITES)
print CANCER_90_METABOLITES
print sorted(CANCER_90_METABOLITES, reverse=True)

CANCER_90_METABOLITE_DEGREE=[]
for met in CANCER_90_METABOLITES:
    CANCER_90_METABOLITE_DEGREE.append(METABO_NT.out_degree(met))

CANCER_90_METABOLITE_DEGREE=filter(lambda x: x<50 and x>35, CANCER_90_METABOLITE_DEGREE)
plt.hist(CANCER_90_METABOLITE_DEGREE)
plt.show()
