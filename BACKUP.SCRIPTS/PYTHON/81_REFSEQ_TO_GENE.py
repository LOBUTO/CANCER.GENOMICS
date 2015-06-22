#OBJECT: RefSeq Protein -> Gene # UNUSABLE AT THE END    
#Based on isoformsFasta object from Dario
#8/5/13

import subprocess
import pickle

DICT_REFSEQ_TO_GENE={}
DICT_GENE_TO_REFSEQ={}
FASTA_FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/isoformsFasta").splitlines()

for file in FASTA_FILES:
    FILE_IN=open("DATABASES/CANCER_DATA/isoformsFasta/%s"%file)
    REF_SEQ_ID=FILE_IN.read().splitlines()[0].split(">")[1].strip()
    DICT_REFSEQ_TO_GENE[REF_SEQ_ID]=file.split(".fasta")[0]
    DICT_GENE_TO_REFSEQ[file.split(".fasta")[0]]=REF_SEQ_ID
    

print DICT_REFSEQ_TO_GENE
print DICT_GENE_TO_REFSEQ
print len(DICT_REFSEQ_TO_GENE)
print len(DICT_GENE_TO_REFSEQ)
print len(FASTA_FILES)

#Save to OBJECT DICT
PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_CANCER_REFSEQ_TO_GENE", "w")
pickle.dump(DICT_REFSEQ_TO_GENE, PICKLE_OUT1)
PICKLE_OUT2=open("DATABASES/OBJECTS/DICT_CANCER_GENE_TO_REFSEQ", "w")
pickle.dump(DICT_GENE_TO_REFSEQ, PICKLE_OUT2)
PICKLE_OUT1.close()
PICKLE_OUT2.close()

#Print all refseq to file
CANCER_REFSEQ=[X for X in DICT_REFSEQ_TO_GENE.iterkeys()]
FILE_OUT=open("DATABASES/CANCER_DATA/ALL_REFSEQ", "w")
for record in CANCER_REFSEQ:
    FILE_OUT.write(record+"\n")
FILE_OUT.close()