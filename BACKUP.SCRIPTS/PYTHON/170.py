#FUNCTION.REMOVE.OVERLAP
#011915
#Removes overlaps of exons in exons.file such 110914.EXON.COORDIANTES



def remove_overlap(ranges): #TEST!!
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

def command_run(hugo):
    
    #######EXON FILE PROCESSING########        
    
    #Filter for hugo of interest
    HUGO=filter(lambda x: x[0]==hugo, EXONS)
    RECORD=list(set([tuple(X[0:6]) for X in HUGO]))
    gene, chrm, strand, start, end, feat=RECORD[0]
    
    #Organize into unique set of tuple ranges
    HUGO_EXONS=list(set([tuple([int(x[6]),int(x[7])]) for x in HUGO]))
    
    #Fix overlaps and sort
    HUGO_EXONS=remove_overlap(HUGO_EXONS)
    
    #######WRITE TO FILE##############
    for record in HUGO_EXONS:
        FILE_OUT.write("\n"+gene+"\t"+chrm +"\t"+strand+"\t"+start+"\t"+end+"\t"+feat+"\t"+str(record[0])+"\t"+str(record[1]))
    

if __name__ == "__main__":
    
    ####Open Exon coordinates file
    FILE_IN=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES")
    EXONS=[X.split("\t") for X in FILE_IN.read().splitlines()[1:]]
    FILE_IN.close()
    
    ####Open File to write
    FILE_OUT=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/011915.EXON.COORDIANTES.FIXED.OVERLAPS","w")
    for i in ["Hugo_Symbol","Chrom","Strand","START","END","FEATURE","FEAT_START"]:
        FILE_OUT.write(i + "\t")
    FILE_OUT.write("FEAT_END")
    
    ####Obtain all Hugos
    HUGOS=list(set([Y[0] for Y in EXONS]))
    
    ####Execute
    for hugo in HUGOS:
        command_run(hugo)
    
    ####Close
    FILE_OUT.close()
    
