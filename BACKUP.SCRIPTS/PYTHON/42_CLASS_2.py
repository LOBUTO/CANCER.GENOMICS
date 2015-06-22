#SCRIPT 2 - FILTERING CLUSTERS

import subprocess

def MAKE_DICT_LIST(FILE_IN): #DICTIONARY - MODIFIED 
    FILE=open(FILE_IN)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split()[0:1]+(i.split()[1]+"|"+i.split()[2]).split()) #MODIFIED
    
    print FILEA
    
    MOL_DICT={} #VERY FAST - DICTIONARY
    for key,val in FILEA:
        MOL_DICT.setdefault(key,[]).append(val)
    
    FILE.close()
    return MOL_DICT

SEP=[]

#STEP 1 - Filter for uniqueness within file, removes same ID within subgroups and across subgroups withing same group
CLUSTERS=subprocess.check_output("ls",cwd="CLASS_PROJECT/GROUPS").splitlines()

for i in CLUSTERS:    
    subprocess.call("sort CLASS_PROJECT/GROUPS/%s | uniq |sort -k2.1 |uniq -u -f1 | sort > temp"%i, shell=True)
    subprocess.call("rm CLASS_PROJECT/GROUPS/%s"%i, shell=True)
    subprocess.call("mv temp CLASS_PROJECT/GROUPS/%s"%i, shell=True)

#STEP 2 - Filter for groups that don't have at least 4 IDs in both subgroups #FOUND 17 CLUSTERS
for i in CLUSTERS:
    CLUSTER_DICT=MAKE_DICT_LIST("CLASS_PROJECT/GROUPS/%s"%i)
    
    #Filter first for cluster that have two subgroups
    if len(CLUSTER_DICT)==2:
        CLUSTER_KEYS=[]
        for j in CLUSTER_DICT.iterkeys():
            CLUSTER_KEYS.append(j)
         
        #Filter for subgroups bigger than 3    
        if len(CLUSTER_DICT[CLUSTER_KEYS[0]])>3 and len(CLUSTER_DICT[CLUSTER_KEYS[1]])>3:
            subprocess.call("cp CLASS_PROJECT/GROUPS/%s CLASS_PROJECT/GROUPS_FILTER"%i, shell=True)
