#miRNA project 4
#Get sequence for mapped mature mirna
#031614

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

FILE_IN1=open("DATABASES/CANCER_DATA/miRNA/031614_mature_coor")
mirna=[x.split("\t") for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/031614_mature_coor_seq","w")
FILE_OUT1.write("Chromosome"+"\t"+"Start"+"\t"+"End"+"\t"+"Strand"+"\t"+"Accession"+"\t"+"Name"+"\t"+"Sequence")
x=0

for record in mirna:
    print record, len(mirna), x
    FILE_OUT1.write("\n"+record[0]+"\t"+record[1]+"\t"+record[2]+"\t"+record[3]+"\t"+record[4]+"\t"+record[5])
    
    #Load chromosome
    chromosome=record[0].strip("chr")
    chromosome_file=open("SOFTWARE/pipeline/data/hsRef/hs_ref_GRCh37.p10_chr%s.fa"%chromosome)
    CHR=chromosome_file.read().splitlines()[1:]
    chromosome_file.close()
    CHR_SEQ="".join(CHR)
    
    #Map - Some are in the negative sense strand
    Strand=record[3]
    if Strand=="+":
        sequence=CHR_SEQ[int(record[1])-1:int(record[2])]
    else:
        sequence=Seq(CHR_SEQ[int(record[1])-1:int(record[2])], generic_dna).reverse_complement()
    
    FILE_OUT1.write("\t"+str(sequence))
    x=x+1
FILE_OUT1.close()