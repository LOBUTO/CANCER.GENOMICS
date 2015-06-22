#miRNA project 5
#Map mutations from BRCA.maf to pri-mirna
#031714

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

FILE_IN1=open("DATABASES/CANCER_DATA/miRNA/031714_BRCA.maf.mirna")
HEADER=FILE_IN1.read().splitlines()[0].split("\t")
FILE_IN1.close()

FILE_IN1=open("DATABASES/CANCER_DATA/miRNA/031714_BRCA.maf.mirna")
BRCA_MAF=[X.split("\t") for X in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/031714_BRCA.maf.mirna.mutations","w")
for head in HEADER:
    FILE_OUT1.write(head+"\t")
FILE_OUT1.write("Mut.Seq")
    
for record in BRCA_MAF:
    STRAND=record[10]
    SEQUENCE=record[12]
    MUTATION=record[5]
    REFERENCE=record[4]
    
    CUT_START=int(record[2])-int(record[8])
    CUT_END=int(record[3])-int(record[8])+1
    
    if STRAND=="+":
        if REFERENCE!="-": # To account for insertions
            MUTATED_SEQ=SEQUENCE[:CUT_START]+MUTATION.lower()+SEQUENCE[CUT_END:]
        else:
            MUTATED_SEQ=SEQUENCE[:CUT_START] + MUTATION.lower() + SEQUENCE[CUT_START:]
        
    else:
        if REFERENCE!="-": # To account for insertions
            SEQUENCE=str(Seq(SEQUENCE, generic_dna).reverse_complement()) #To account for negative sense mirna
            MUTATED_SEQ=SEQUENCE[:CUT_START]+MUTATION.lower()+SEQUENCE[CUT_END:] #Place mutation with respect to positive strand
            MUTATED_SEQ=str(Seq(MUTATED_SEQ, generic_dna).reverse_complement())
        else:
            SEQUENCE=str(Seq(SEQUENCE, generic_dna).reverse_complement()) #To account for negative sense mirna
            MUTATED_SEQ=SEQUENCE[:CUT_START]+MUTATION.lower()+SEQUENCE[CUT_START:] #Place mutation with respect to positive strand
            MUTATED_SEQ=str(Seq(MUTATED_SEQ, generic_dna).reverse_complement())
    
    print record, MUTATED_SEQ
    #Write out
    FILE_OUT1.write("\n")
    for info in record:
        FILE_OUT1.write(info+'\t')
    FILE_OUT1.write(MUTATED_SEQ)

FILE_OUT1.close()
        
    