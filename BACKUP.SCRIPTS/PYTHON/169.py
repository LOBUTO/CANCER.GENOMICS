#CHEN.WIGFIX PROCESSING
#120814
#Processing of chen replication time file (from "Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes")
#This step preps the file for annovar annoation (hg18)
#STEP 1 of CHEN.REPLICATION.GENE.ASSIGN.sh

import sys

def command_run (chen_wigfix, file_out):
    
    with open(chen_wigfix) as infile:
        
        #Skip header lines
        for _ in xrange(2):
            next(infile)
        
        #Process 
        for line in infile:
            
            #For every "fixed line" update step change
            if "fixed" in line:
                line_split=line.split()
                chrom=line_split[1].split("=chr")[1]
                start=int(line_split[2].split("=")[1])
                step=int(line_split[3].split("=")[1])
                
            #Use step change to write replication times per reference position and write to file
            else:
                for record in [chrom, str(start), str(start), "0","0"]:
                    file_out.write(record+"\t")
                file_out.write(line)
                start=start+step
                
if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT=sys.argv[1], open(sys.argv[2],"w")
    #FILE_IN1="DATABASES/CANCER_DATA/REP.TIME/S50.100kbp.Non-overlappingWindows.hg18.wig"
    #FILE_OUT=open("LOGS/temp.1","w")
    
    ###Process
    command_run(FILE_IN1, FILE_OUT)
    FILE_OUT.close()