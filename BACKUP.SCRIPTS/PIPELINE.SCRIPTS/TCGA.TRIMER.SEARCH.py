#####TCGA.TRIMER.SEARCH.py######
#123014
#Find TRIMERS for nt positions present in TCGA.maf pre-processed files by Rscript
#Done with human genome build 37

import sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import product

def melt_temp_dict():
    #Create melting temperature dictionary for all trimers
    
    #Create trimer list
    all_trimers=list(product(["A","C","G","T","N"], repeat=3))
    all_trimers=["".join(x) for x in all_trimers]
    ACGT_trimers=filter(lambda x: "N" not in x, all_trimers)
    N_trimers=filter(lambda x: "N" in x, all_trimers)
    
    #Obtain melting temperature
    trimer_dict=dict((x,mt.Tm_Wallace(Seq(x))) for x in ACGT_trimers)
    
    #Add "NA" for "NNN" temperatures
    n_trimer_dict=dict((x,"NA") for x in N_trimers)
    trimer_dict.update(n_trimer_dict)
    
    #Return
    return trimer_dict

def fasta_dict(fasta_file):
    #Create the fasta dictionary with chromosome entries as keys and nt strings as values
    fasta_dict=SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return fasta_dict

def command_run(tcga, FILE_OUT, chrom_dict, melt_dict):
    
    #Open file out
    file_out=open(FILE_OUT, "w")
    
    #Process tcga file
    with open(tcga) as infile:
        
        #Write output header
        header=infile.readline().split("\t")
        for line in header[:-1]:
            file_out.write(line + "\t")
        file_out.write("TRIMER"+"\t"+"MT"+"\t")
        file_out.write(header[-1])
        
        #Process and assign trimer
        for line in infile:
            
            record=line.split("\t")
            
            #Filter out non-single substituions
            if len(record[6])==1:
                
                #Assign trimer
                chrom, pos=record[3], record[0]
                trimer=str(chrom_dict[chrom].seq[int(pos)-2:int(pos)+1])
                mt=str(melt_dict[trimer])
                
                #Write to file
                for element in record[:-1]:
                    file_out.write(element + "\t")
                file_out.write(trimer + "\t" + mt + "\t")
                file_out.write(record[-1])
    
    #Close file out
    file_out.close()
            
if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT, HUMAN_FASTA=sys.argv[1], sys.argv[2], sys.argv[3]
    print "Done loading files"
    
    ####Load fasta file as dict
    FASTA_FILE=open(HUMAN_FASTA) #As in "DATABASES/CANCER_DATA/1000GENOME/2013/human_g1k_v37.fasta"
    chrom_dict=fasta_dict(FASTA_FILE)
    FASTA_FILE.close()
    print "Done constructing fasta dictionary"
    
    ####Create trimer melting temperature dictionary
    melt_dict=melt_temp_dict()
    
    ###Process
    print "Running trimer search on maf file"
    command_run(FILE_IN1, FILE_OUT, chrom_dict, melt_dict)