#DESPRAT.WIGFIX PROCESSING
#011215
#Processing of desprat replication files 

import sys

def command_run (chen_wigfix, file_out):
    
    #Loop to obtain rep times and store in dict
    with open(chen_wigfix) as infile:
        
        #Skip header lines
        for _ in xrange(2):
            next(infile)
        
        #Process 
        for line in infile:
            
            #Obtain attributes
            line_split=line.split("\t")
            chr=line_split[0].split("chr")[1]
            pos=line_split[1]
            rt=line_split[3]
            
            #Write to file
            for i in [chr, pos, pos, "0","0",rt]:
                file_out.write(i+"\t")

if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT=sys.argv[1], open(sys.argv[2],"w")
    
    ###Process
    command_run(FILE_IN1, FILE_OUT)
    FILE_OUT.close()