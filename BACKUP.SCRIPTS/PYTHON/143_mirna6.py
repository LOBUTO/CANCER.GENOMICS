#miRNA project 6
#Map mutations from BRCA.maf to mature miRs
#031714
#NOTE: Will only map if mutation falls completely within mature sequence, otherwise will discard
#Will write function for it

def mirna_maf_mut(maf_mirna, out_file): #Will take maf.mirna or maf.mature as long as is in the 031714_BRCA.maf.mirna format
    COUNT=0 #To keep track
    
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    
    FILE_IN1=open(maf_mirna)
    HEADER=FILE_IN1.read().splitlines()[0].split("\t")
    FILE_IN1.close()
    
    FILE_IN1=open(maf_mirna)
    MAF=[X.split("\t") for X in FILE_IN1.read().splitlines()[1:]]
    FILE_IN1.close()
    
    FILE_OUT1=open(out_file,"w")
    for head in HEADER:
        FILE_OUT1.write(head+"\t")
    FILE_OUT1.write("Mut.Seq")
        
    for record in MAF:
        
        #First check if mutations is within boudaries of mirna (this is mainly a check for maf.mature files)
        CHECK_UPSTREAM=int(record[2])-int(record[8])
        CHECK_DOWNSTREAM=int(record[9])-int(record[3])
        if CHECK_UPSTREAM>=0 and CHECK_DOWNSTREAM>=0:
            COUNT=COUNT+1
            
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
            
            #Write out
            FILE_OUT1.write("\n")
            for info in record:
                FILE_OUT1.write(info+'\t')
            FILE_OUT1.write(MUTATED_SEQ)
    
    print maf_mirna.split("_")[-1].split(".")[0].strip("maf"),"\t",maf_mirna.split("_")[-1].split(".")[1],"\t",COUNT
    
    FILE_OUT1.close()

SEPARATOR=[]

#mirna_maf_mut("DATABASES/CANCER_DATA/miRNA/031814_BRCA.maf.mature", "DATABASES/CANCER_DATA/miRNA/031814_BRCA.maf.mature.mutations")    

import subprocess

#Process for all cancers that have mirna data
FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/miRNA/").splitlines()    
FILES=filter(lambda x: "maf.mature" in x[-10:] or "maf.mirna" in x[-9:], FILES)

for record in FILES:
    mirna_maf_mut("DATABASES/CANCER_DATA/miRNA/%s"%record,"DATABASES/CANCER_DATA/miRNA/%s.mutations"%record)

