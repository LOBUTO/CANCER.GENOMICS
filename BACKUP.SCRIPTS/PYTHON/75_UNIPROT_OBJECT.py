#MAKE UNIPROT - GENE DICTIONARY OBJECT
import pickle

#FOR ALL UNIPROTS IN DATABASE
FILE_IN1=open("/Users/jzamalloa/Desktop/FOLDER/LOGS/uniprot_to_gene.text").read().splitlines()
print FILE_IN1[:3]
DICT_UNIPROT_TO_GENES={}
for record in FILE_IN1:
    UNIPROT=record.split("\t")[0].strip()
    GENES=record.split("\t")[1].strip().split()
    DICT_UNIPROT_TO_GENES[UNIPROT]=GENES

PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_GENES.pi","w")
pickle.dump(DICT_UNIPROT_TO_GENES,PICKLE_OUT1)
print DICT_UNIPROT_TO_GENES["Q96QU6"]
print DICT_UNIPROT_TO_GENES["P30493"]
print DICT_UNIPROT_TO_GENES["P33765"]  
print DICT_UNIPROT_TO_GENES["A2A322"]   
print DICT_UNIPROT_TO_GENES["Q8IZR3"]

#FOR ONLY REVIEWED UNIPROTS IN DATABASE
FILE_IN2=open("/Users/jzamalloa/Desktop/FOLDER/LOGS/uniprot_to_gene_reviewed.text")
DICT_UNIPROT_TO_GENES_REVIEWED={}
for record in FILE_IN2:
    UNIPROT=record.split("\t")[0].strip()
    GENES=record.split("\t")[1].strip().split()
    DICT_UNIPROT_TO_GENES_REVIEWED[UNIPROT]=GENES

PICKLE_OUT2=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_GENES_REVIEWED.pi", "w")
pickle.dump(DICT_UNIPROT_TO_GENES_REVIEWED, PICKLE_OUT2)