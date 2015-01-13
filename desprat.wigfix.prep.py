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
            rt=line_split[3].rstrip("\r\n")
            
            #Write to file
            for i in [chr, pos, pos, "0","0",rt]:
                file_out.write(i+"\t")
            file_out.write("\n")    
            
if __name__ == "__main__":
    
    ####Load files
    FILE_IN1, FILE_OUT=sys.argv[1], open(sys.argv[2],"w")
    #FILE_IN1="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/REP.TIME/Supplementary_Data_090109/25MS_wig.txt"
    #FILE_OUT=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/LOGS/temp.1","w")
    
    ###Process
    command_run(FILE_IN1, FILE_OUT)
    FILE_OUT.close()
    