#Processing PHASTCONS without QUEUE (PARALLEL, no loading files)

#Need to assign process wig files into order to obtain CDS positions only information (Processing of Phastcons files)
#Need wigfix files and exon.coordinates files for this

import multiprocessing as mp
import subprocess, time, threading

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

def foo(item):
        print item, time.ctime()
        threading.Timer(5, foo, [item]).start()

def command_run(CHROMOSOME):
    
    
    print "chromosome", CHROMOSOME
    
    #######EXON FILE PROCESSING########        
    
    #Filter for chromosome of interest
    CHRM_EXONS=filter(lambda x: x[1]==CHROMOSOME, EXONS)
    
    #Organize into unique set of tuple ranges
    CHRM_EXONS=list(set([tuple([int(x[3]),int(x[4])]) for x in CHRM_EXONS]))
    
    #Fix overlaps and sort
    CHRM_EXONS=remove_overlap(CHRM_EXONS)
    
    print "Processed EXONS"
    #######################    
    
    #######OPEN CHROMOSOME SPECIFIC FILE TO WRITE TO
    FILE_OUT1=open("DATABASES/PHASTCONS/45.PLACENTAL.PROCESSED/120414.CHRM."+CHROMOSOME+".FILTERED.EXONS.45.PARALLEL", "w")
    FILE_OUT1.write("Chrom"+"\t"+"Position"+"\t"+"PHAST"+"\n")
    
    #######LOAD CURRENT PHASTCON######
    
    #Get length of file (WIG FILE)
    FILE="DATABASES/PHASTCONS/45.PLACENTAL/chr"+ CHROMOSOME +".phastCons46way.placental.wigFix"
    file_length=int(subprocess.check_output(["wc", "-l", FILE]).split()[0])
    
    #Print out count every 10 seconds
    count=1
    foo([CHROMOSOME,float(count)/file_length])
    
    with open(FILE) as infile:
        for line in infile:
            
            #For every fixed step break line we update the START and STEP record
            if "fixedStep" in line:
                START=int(line.split()[2].split("=")[1])
                STEP=int(line.split()[3].split("=")[1])
                COUNT=0
                
            #If no fixed step found in line then we continue with steps of current fixed step
            else:
                
                #Check that position (START+COUNT) exists in chromosomal exon range for exons, if it is then write it
                if binarySearch(CHRM_EXONS, START+COUNT):
                    
                    FILE_OUT1.write(CHROMOSOME +"\t"+str(START+COUNT)+"\t"+line)
                
                #Update count - That is, updating nucleotide step we are in right now
                COUNT=COUNT+STEP
                
            count=count+1
    
    ########CLOSE FILE
    FILE_OUT1.close()
    
#####

if __name__ == "__main__":
    
    ####Open Exon coordinates file
    FILE_IN=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES")
    EXONS=[X.split("\t") for X in FILE_IN.read().splitlines()[1:]]
    FILE_IN.close()
    
    ####Set up chromosomes
    chromosomes=[str(x) for x in range(1,23)] + ["X","Y"]
    
    ######Parallelize job
    pool = mp.Pool()
    pool.map(command_run, chromosomes)
    
    #When done stop
    pool.close()
    pool.join()
    