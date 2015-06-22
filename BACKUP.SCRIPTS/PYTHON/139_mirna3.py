#miRNA project 3
#Map sequence to primirna using coordinates from mirbase
#031614

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#Map
FILE_IN1=open("DATABASES/CANCER_DATA/miRNA/031614_primirna_coor")
mirna=[x.split("\t") for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/031614_primirna_coor_seq","w")
FILE_OUT1.write("Synonyms"+"\t"+"Symbol"+"\t"+"Chromosome"+"\t"+"Start"+"\t"+"End"+"\t"+"Strand"+"\t"+"Accession"+"\t"+"Sequence")
x=0
for record in mirna:
    print record, len(mirna), x
    FILE_OUT1.write("\n"+record[0]+"\t"+record[1]+"\t"+record[2]+"\t"+record[3]+"\t"+record[4]+"\t"+record[5]+"\t"+record[6])
    
    #Load chromosome
    chromosome=record[2].strip("chr")
    chromosome_file=open("SOFTWARE/pipeline/data/hsRef/hs_ref_GRCh37.p10_chr%s.fa"%chromosome)
    CHR=chromosome_file.read().splitlines()[1:]
    chromosome_file.close()
    CHR_SEQ="".join(CHR)
    
    #Map - Some are in the negative sense strand
    Strand=record[5]
    if Strand=="+":
        sequence=CHR_SEQ[int(record[3])-1:int(record[4])]
    else:
        sequence=Seq(CHR_SEQ[int(record[3])-1:int(record[4])], generic_dna).reverse_complement()
    
    FILE_OUT1.write("\t"+str(sequence))
    x=x+1
FILE_OUT1.close()