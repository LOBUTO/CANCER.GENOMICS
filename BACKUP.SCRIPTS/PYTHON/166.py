#MAP DATA FROM ONE THOUSAND GENOME TO TRI-MERS - DONE!! - NEED TO FIX SINCE TRIMERS ARE NOT ACCURATE DUE TO LACK OF COMPLEMENT CONSIDERATION!!!!
#121914
#IMPLEMENTED:
#    Obtaining melting temperature per primer 122014
#    Classifying each position as exon or intron

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

def remove_overlap(ranges):
    result = []
    current_start = -1
    current_stop = -1 
    
    #Sort by first item
    for start, stop in sorted(ranges):
        
        if start > current_stop:
            # this segment starts after the last segment stops
            # just add a new segment
            result.append( (start, stop) )
            current_start, current_stop = start, stop
        
        else:
            
            #Consider max stop
            current_stop = max(current_stop, stop)
            
            # segments overlap, replace
            result[-1] = (current_start, current_stop)

    return result
    
def binarySearch(alist, item):
    if len(alist) == 0:
        return False
    
    else:
        midpoint = len(alist)//2
        if alist[midpoint][0]<=item<=alist[midpoint][1]:
            return True
    
        else:
            if item<alist[midpoint][0]:
                return binarySearch(alist[:midpoint],item)
            else:
                return binarySearch(alist[midpoint+1:],item)

def exon_dict(EXON_LIST):
    #Create exon dictionary per chromosome
    
    #Set up chromosomes
    chromosomes=[str(x) for x in range(1,23)] + ["X","Y"]
    
    #Create storage dictionary
    EXON_DICT=dict()
    
    #######EXON FILE PROCESSING########        
    for chrom in chromosomes:
        
    
        #Filter for chromosome of interest
        CHRM_EXONS=filter(lambda x: x[1]==chrom, EXONS)
        
        #Organize EXONS into unique set of tuple ranges (FEAT_START and FEAT_END)
        CHRM_EXONS=list(set([tuple([int(x[6]),int(x[7])]) for x in CHRM_EXONS]))
        
        #Fix overlaps and sort
        CHRM_EXONS=remove_overlap(CHRM_EXONS)
        
        #Store filtered chromosomal coordinates into dictionary
        EXON_DICT[chrom]=CHRM_EXONS
        
    ####################################
    
    print "Processed EXONS"
    
    #Return
    return EXON_DICT

def command_run(thousand_file, chrom_dict, file_out, melt_dict, exon_dict):
    
    with open(thousand_file) as infile:
        for line in infile:
            
            #Filter out non single substitutions
            record=line.split("\t")
            if (len(record[3])==1 and len(record[4])==1 and record[2]=="SNP"):
                
                #Parse
                chrom, pos, ref, alt, maf=record[0:2]+record[3:5]+[record[10]]
                trimer=str(chrom_dict[chrom].seq[int(pos)-2:int(pos)+1])
                mt=str(melt_dict[trimer])
                region=str(binarySearch(exon_dict[chrom], int(pos)))
                
                #Write out to file
                for i in [chrom, pos, ref, alt, trimer, mt, region]:
                    file_out.write(i +"\t")
                file_out.write("\t"+maf)         

if __name__ == "__main__":
    
    ####Load fasta file as dict
    FASTA_FILE=open("DATABASES/CANCER_DATA/1000GENOME/2013/human_g1k_v37.fasta")
    chrom_dict=fasta_dict(FASTA_FILE)
    FASTA_FILE.close()

    ####Create trimer melting temperature dictionary
    melt_dict=melt_temp_dict()
    
    ####Load Exon coordinates file
    FILE_IN=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES")
    EXONS=[X.split("\t") for X in FILE_IN.read().splitlines()[1:]]
    FILE_IN.close()
    
    ####Convert Exon coordiates into dict per chromosome
    EXON_DICT=exon_dict(EXONS)
    
    ####Open 1000G processed file
    THOUSAND_FILE="PIPELINES/METABOLIC.DRIVERS/TABLES/120814.THOUSAND.SNP.ALL.CNV"
    
    ####Open file to write to
    FILE_OUT=open("PIPELINES/METABOLIC.DRIVERS/TABLES/122014.THOUSAND.SNP.TRIMERS.PLUS", "w")
    FILE_OUT.write("Chrom" +"\t" + "Position" +"\t" + "REF" + "\t" + "ALT" + "\t" + "TRIMER" + "\t" + "MT" +"\t" +"REGION"+"\t"+
                    "AF" +"\n")
    
    ####Fire main program
    command_run(THOUSAND_FILE, chrom_dict, FILE_OUT, melt_dict, EXON_DICT)
    FILE_OUT.close()

