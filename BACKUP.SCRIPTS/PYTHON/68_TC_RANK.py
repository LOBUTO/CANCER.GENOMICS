######TANIMOTO COEFFICIENT CALCULATION###### - Wrapper for babel
####Calculations are based on comparisson of a file containing a single line of smiles vs a a file containing either a list of smiles or
####.sdf molecules (MAKE SURE THEY ARE WELL FORMATTED). The output will produce a list ranked from highest to lowest score
####This script depends on a machine with OpenBabel installed (http://openbabel.org/wiki/Python), IMPORTANT make sure path is correct in script
####IMPORTANT - Default fingerprint type is FP2, could choose instead: FP3, FP4 or MACCS
####6/13/13 - JOSE A. ZAMALLOA

def TC_RANK(QUERY, LIST, OUTPUT, FP="FP2"): #QUERY = Single smiles code(no name),
                                            #LIST = To compare against, in the format "SMILES    NAME" per line
                                            #OUTPUT file
    import subprocess
    
    TERMINAL_RESULT=subprocess.Popen(["-c", "/usr/local/bin/babel %s %s -ofpt -xf%s"%(QUERY, LIST, FP)], stdout=subprocess.PIPE, shell=True)
    TERMINAL_RESULT=TERMINAL_RESULT.communicate()[0].splitlines()

    #Filter out troublesome lines
    TERMINAL_RESULT=filter(lambda x:len(x)>3 and x!="Possible superstructure of first mol", TERMINAL_RESULT)

    #Rank in descending order
    TC_RANK=[]
    for tc in TERMINAL_RESULT:
        TC_RANK.append([float(tc.split("=")[-1]), tc[0:(len(tc)-len(tc.split("=")[-1]))]])
    
    TC_RANK=sorted(TC_RANK, reverse=True)    
    
    #Write file
    FILE_OUT=open(OUTPUT,"w")
    for line in TC_RANK:
        FILE_OUT.write(line[1]+" "+str(line[0])+"\n")
    FILE_OUT.close()
    

#####TESTING#######
TEST_CASE=["/Users/jzamalloa/Desktop/FOLDER/LOGS/HMDB10366.smi", 
        "/Users/jzamalloa/Desktop/FOLDER/LOGS/END_SMILES.smi",
        "/Users/jzamalloa/Desktop/FOLDER/LOGS/OUT"]

TC_RANK(TEST_CASE[0], TEST_CASE[1], TEST_CASE[2])