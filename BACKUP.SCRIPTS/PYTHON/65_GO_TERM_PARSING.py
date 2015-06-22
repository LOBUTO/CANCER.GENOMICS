#PROCESS GO_TERM FINDER FILES AND OBTAIN GO TERMS FOUND IN EACH FILE

import subprocess
import os

def GO_TERM_FINDER_PARSING(FILE_IN): #Removes non-table lines of GO_TERM_FINDER output files, make sure FILE_IN has full path with respect
                                    #to working directory
    subprocess.call('sed "1,11d" %s > TEST'%FILE_IN, shell=True)
    subprocess.call("rm %s"%FILE_IN, shell=True)
    subprocess.call("mv TEST %s"%FILE_IN, shell=True)

def GO_TERM_FINDER_PARSED_TERMS(FILE_IN, FILE_OUT): #Takes files from GO_TERM_FINDER_PARSING and produces FILE_OUT with single_spaced list 
                                                    #of GO TERMS
                                                    #IMPORTANT, FILE_IN needs full path with respect to working directory
    FILE_UNREAD=open(FILE_IN)
    FILE_READ=FILE_UNREAD.read().splitlines()
    FILE_UNREAD.close()
    
    OUTPUT=open(FILE_OUT, "w")
    for line in FILE_READ[1:]:
        OUTPUT.write(line.split("\t")[0]+"\n")
    OUTPUT.close()

"""    
GO_FILES=subprocess.check_output("ls", cwd="NETWORK/GO_ENRICHMENT").splitlines()

for file in GO_FILES:
    file="NETWORK/GO_ENRICHMENT/"+file
    GO_TERM_FINDER_PARSING(file)
    GO_TERM_FINDER_PARSED_TERMS(file, file+".TERMS")
"""    

#Count number of terms (lines) in each .GO.TERMS file
os.chdir("NETWORK/GO_ENRICHMENT")
FOR_FILES=subprocess.Popen(["-c", "ls *.TERMS"], stdout=subprocess.PIPE, shell=True)
FOR_FILES=FOR_FILES.communicate()[0].splitlines()
print FOR_FILES

for file in FOR_FILES:
    print subprocess.check_output(["wc", "-l", "%s"%file]).splitlines()
