#####TCGA.TRINUCLEOTIDE.SEARCH.py###### 
#123014
#Find TRIMERS for nt positions present in TCGA.maf pre-processed files by Rscript
#Done with human genome build 37

import sys
from Bio import SeqIO

def fasta_dict(fasta_file):
    #Create the fasta dictionary with chromosome entries as keys and nt strings as values
    fasta_dict=SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return fasta_dict

def command_run(tcga, file_out, chrom_dict):
    
    with open(tcga) as infile:
        
        #Write output header
        header=infile.readline().split("\t")
        for line in header[:-1]:
            file_out.write(line + "\t")
        file_out.write("TRIMER"+"\t")
        file_out.write(header[-1])
        
        #Process and assign trimer
        for line in infile:
            
            record=line.split("\t")
            
            #Filter out non-single substituions
            if len(record[4])==1:
                
                #Assign trimer
                chrom, pos=record[1], record[2]
                trimer=str(chrom_dict[chrom].seq[int(pos)-2:int(pos)+1])
                
                #Write to file
                for element in record[:-1]:
                    file_out.write(element + "\t")
                file_out.write(trimer + "\t")
                file_out.write(record[-1])
            
if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT, HUMAN_FASTA=sys.argv[1], sys.argv[2], sys.argv[3]
    
    ####Load fasta file as dict
    FASTA_FILE=open(HUMAN_FASTA) #As in "DATABASES/CANCER_DATA/1000GENOME/2013/human_g1k_v37.fasta"
    chrom_dict=fasta_dict(FASTA_FILE)
    FASTA_FILE.close()
    
    ###Process
    command_run(FILE_IN1, FILE_OUT, chrom_dict)