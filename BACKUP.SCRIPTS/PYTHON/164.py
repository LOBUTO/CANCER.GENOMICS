#Processing PHASTCONS with QUEUE (PARALLEL, no loading files)

#Need to assign process wig files into order to obtain CDS positions only information (Processing of Phastcons files)
#Need wigfix files and exon.coordinates files for this

import multiprocessing as mp
import subprocess

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

def load_wigs (CHROMOSOMES):
    INTERNAL_DICT=dict((chrm,[]) for chrm in CHROMOSOMES)
    
    for chrm in CHROMOSOMES:
        print chrm
        FILE_IN=open("DATABASES/PHASTCONS/45.TRACK/chr"+ chrm +".phastCons46way.wigFix")
        COOR=[X for X in FILE_IN.read().splitlines()]
        FILE_IN.close()
        INTERNAL_DICT[chrm]=COOR
        
    return INTERNAL_DICT

def command_run(CHROMOSOME, q):
    print "chromosome", CHROMOSOME
    
    #######EXON FILE PROCESSING########        
    
    #Filter for chromosome of interest
    CHRM_EXONS=filter(lambda x: x[1]==CHROMOSOME, EXONS)
    
    #Organize into unique set of tuple ranges
    CHRM_EXONS=list(set([tuple([int(x[3]),int(x[4])]) for x in CHRM_EXONS]))
    
    #Fix overlaps and sort
    CHRM_EXONS=remove_overlap(CHRM_EXONS)
    
    print "Processed EXONS"
    ##################################    
    
    #######LOAD CURRENT PHASTCON######
    
    #Get length of file
    FILE="DATABASES/PHASTCONS/45.TRACK/chr"+ CHROMOSOME +".phastCons46way.wigFix"
    file_length=int(subprocess.check_output(["wc", "-l", FILE]).split()[0])
    
    count=1
    with open(FILE) as infile:
        for line in infile:
        
            print CHROMOSOME, (float(count)/file_length)
            
            #For every fixed step break line we update the START and STEP record
            if "fixedStep" in line:
                START=int(line.split()[2].split("=")[1])
                STEP=int(line.split()[3].split("=")[1])
                DICT_COUNT={CHROMOSOME:0}
                #COUNT=0
                
            #If no fixed step found in line then we continue with steps of current fixed step
            else:
                
                #Check that position (START+COUNT) exists in chromosomal exon range for exons, if it is then write it
                if binarySearch(CHRM_EXONS, START+DICT_COUNT[CHROMOSOME]):
                    q.put(CHROMOSOME +"\t"+str(START+DICT_COUNT[CHROMOSOME])+"\t"+line)
                
                #Update count - That is, updating nucleotide step we are in right now
                DICT_COUNT[CHROMOSOME]=DICT_COUNT[CHROMOSOME]+STEP
                
            count=count+1

def listener(q):
    '''listens for messages on the q, writes to file. ''' 
    
    FILE_OUT1=open("DATABASES/PHASTCONS/12.02.14.CHRM.FILTERED.EXONS.45.PARALLEL", "w")
    FILE_OUT1.write("Chrom"+"\t"+"Position"+"\t"+"PHAST"+"\n")
    
    while 1:
        m = q.get()
        if m == 'kill':
            break
        FILE_OUT1.write(m)
        FILE_OUT1.flush()
    FILE_OUT1.close()
    
#####

if __name__ == "__main__":
    ####Open Exon coordinates file
    FILE_IN=open("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES")
    EXONS=[X.split("\t") for X in FILE_IN.read().splitlines()[1:]]
    FILE_IN.close()
    
    ####Open WIG files a create dictionary to process
    chromosomes=[str(x) for x in range(1,23)] + ["X","Y"]
    
    ######Parallelize job
    
    #Set up queue
    manager = mp.Manager()
    q = manager.Queue()    
    pool = mp.Pool()
    print "queue set up"
    
    #Start listener
    watcher = pool.apply_async(listener, (q,))
    print "listener ready"
    
    #Fire workers
    jobs = []
    for chrm in chromosomes:
        job=pool.apply_async(command_run, (chrm,q))
        jobs.append(job)
    
    # collect results from the workers through the pool result queue
    for job in jobs: 
        job.get()
    
    #When done kill all
    q.put("kill")
    pool.close()
    pool.join()
    