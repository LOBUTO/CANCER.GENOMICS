#CHEN.WIGFIX PROCESSING
#120814
#Processing of chen replication time file (from "Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes")
#This step preps the file for annovar annoation (hg18)
#STEP 1 of CHEN.REPLICATION.GENE.ASSIGN.sh

import sys
import numpy as np

def mid_rep (mid_vector):
    #Do middle average position of rep time to fill in gaps, logical since we are covering neighboring regions of chromosome that
    # are supposed to have similar replication times
    
    end_vector=[]
    for i in range(len(mid_vector)-1):
        end_vector.append(mid_vector[i])
        end_vector.append(  np.mean([mid_vector[i],mid_vector[i+1]]) )
    end_vector.append(mid_vector[-1])
    
    #Return
    return end_vector
    
def command_run (chen_wigfix, file_out):
    
    #Open storage dictionary
    chrm_dict={}
    
    #Loop to obtain rep times and store in dict
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
                
                #For modified mid vector storage in dict
                chrm_dict[chrom]={"pos":[], "mid_vector":[]}
                
            #Use step change to write replication times per reference position and write to file
            else:
                #Apply mid replication time function
                chrm_dict[chrom]["mid_vector"].append(float(line.rstrip("\n")))
                chrm_dict[chrom]["pos"].append(start)
                start=start+step
                
    #Obtain mid rep times per chromosomal vector - Applied extension twice
    for chrm in chrm_dict.keys():
        chrm_dict[chrm]["mid_vector"]=   mid_rep(mid_rep(mid_rep(chrm_dict[chrm]["mid_vector"]))) #RT are floats
        chrm_dict[chrm]["pos"]=[int(x) for x  in    mid_rep(mid_rep(mid_rep(chrm_dict[chrm]["pos"])))] #Because position are integers
    
    #Write to file
    for chrm in chrm_dict.keys():
        
        for record in xrange(len(chrm_dict[chrm]["pos"])):
            file_out.write(chrm+"\t"+ str(chrm_dict[chrm]["pos"][record])+"\t"+ str(chrm_dict[chrm]["pos"][record])+"\t"+
                           "0"+"\t"+"0"+"\t"+
                           str(chrm_dict[chrm]["mid_vector"][record])+"\n")
                
if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT=sys.argv[1], open(sys.argv[2],"w")
    #FILE_IN1="DATABASES/CANCER_DATA/REP.TIME/S50.100kbp.Non-overlappingWindows.hg18.wig"
    #FILE_OUT=open("LOGS/temp.1","w")
    
    ###Process
    command_run(FILE_IN1, FILE_OUT)
    FILE_OUT.close()