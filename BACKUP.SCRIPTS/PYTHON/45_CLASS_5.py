#SCRIPT 5 - Get SDPAs

import subprocess
import os

def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

os.chdir("CLASS_PROJECT")

ALN=subprocess.check_output("ls", cwd="ALIGNMENT_PARSED").splitlines()

for i in ALN:
    #Get both GO
    subprocess.call('''awk 'BEGIN {FS = "|"} ; {print $1}' ALIGNMENT_PARSED/%s |sort |uniq > LOGS/temp'''%i, shell=True)
    GO=open("LOGS/temp")
    GO_TERMS=[]
    for j in GO: GO_TERMS.append(j[:-1])
    GO.close()
    GO1=GO_TERMS[0]
    GO2=GO_TERMS[1]
    
    #Get SDPAs
    CLUSTER_IN=open("ALIGNMENT_PARSED/%s"%i)
    CLUSTER=[]
    for l in CLUSTER_IN: CLUSTER.append(l.split())
    CLUSTER_IN.close()
    
    FILE_OUT=open("SDPA/%s.sdpA"%i, "w")
    SEQ_LEN=len(CLUSTER[0][1])
    X=0
    while X!=SEQ_LEN:
        SDPA_A=[]
        SDPA_B=[]
        for m in CLUSTER:
            if m[0].split("|")[0]==GO1:
                SEQ_A=m[1]
                if SEQ_A[X] not in SDPA_A:
                    SDPA_A.append(SEQ_A[X])
            elif m[0].split("|")[0]==GO2:
                SEQ_B=m[1]
                if SEQ_B[X] not in SDPA_B:
                    SDPA_B.append(SEQ_B[X])
        
        SDPA_A=unique(SDPA_A)
        SDPA_B=unique(SDPA_B)
        
        if len(SDPA_A)==1==len(SDPA_B) and SDPA_A[0]!=SDPA_B[0] and SDPA_A[0]!="-" and SDPA_B[0]!="-":
            FILE_OUT.write(str(X)+"\n")
        X=X+1
    
    FILE_OUT.close()