#Processing PHASTCONS, no parallel unclean (loads file every time)

#Need to assign process wig files into order to obtain CDS positions only information (Processing of Phastcons files)
#Need wigfix files and exon.coordinates files for this

#Try with first wigfix file
#Break file into chunks depending on fixedstep start

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

####Open file to write to
FILE_OUT1=open("DATABASES/PHASTCONS/12.02.14.CHRM.FILTERED.EXONS.45", "w")
FILE_OUT1.write("Chrom"+"\t"+"Position"+"\t"+"PHAST")

####Loop through chromosomes
CHROMOSOMES=[str(x) for x in range(1,23)] + ["X","Y"]

for chrm in CHROMOSOMES:
    print chrm
    
    FILE_IN1=open("DATABASES/PHASTCONS/45.TRACK/chr"+ chrm +".phastCons46way.wigFix")
    COOR=[X for X in FILE_IN1.read().splitlines()]
    FILE_IN1.close()
    
    count=1

    #######EXON FILE########
    #Open Exon coordinates file
    FILE_IN=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES")
    EXONS=[X.split("\t") for X in FILE_IN.read().splitlines()[1:]]
    FILE_IN.close()
    
    #Filter for chromosome of interest
    EXONS=filter(lambda x: x[1]==chrm, EXONS)
    
    #Organize into unique set of tuple ranges
    EXONS=list(set([tuple([int(x[3]),int(x[4])]) for x in EXONS]))
    
    #Fix overlaps and sort
    EXONS=remove_overlap(EXONS)
    
    print "Processed EXONS"
    #######################
    
    for record in COOR:
    
        print chrm, (float(count)/len(COOR))
        
        #For every fixed step break line we update the START and STEP record
        if "fixedStep" in record:
            START=int(record.split()[2].split("=")[1])
            STEP=int(record.split()[3].split("=")[1])
            COUNT=0
            
        #If no fixed step found in line then we continue with steps of current fixed step
        else:
            
            #Check that position (START+COUNT) exists in chromosomal exon range for exons, if it is then write it
            if binarySearch(EXONS, START+COUNT):
                FILE_OUT1.write("\n"+ chrm +"\t"+str(START+COUNT)+"\t"+record)
            
            #Update count - That is, updating nucleotide step we are in right now
            COUNT=COUNT+STEP
            
        count=count+1

####Close file
FILE_OUT1.close()