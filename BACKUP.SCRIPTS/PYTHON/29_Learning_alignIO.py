from Bio import AlignIO

FILE=open("capra_bioinf08_data/1.1.1.1_group0_PF00106_30-85~1.1.1.35_group0_PF00106_30-85.aln")
OUTPUT=open("BIOPYTHON/TEST1","w")

alignments=AlignIO.parse(FILE,"clustal")
AlignIO.write(alignments, OUTPUT,"stockholm") #stockholm format presents them as lines

FILE.close()
OUTPUT.close()
