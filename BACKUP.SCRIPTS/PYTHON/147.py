import subprocess
import numpy as np

FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS").splitlines()
FILES=filter(lambda x: "NOMEGA.T.0.5_5C" in x, FILES)
print FILES

for record in FILES:
    FILE_IN=open("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS"+"/"+record)
    KEGGS=[y.split() for y in FILE_IN.read().splitlines()]
    FILE_IN.close()
    
    ALL_KEGGS=[]
    AVERAGE_CLUSTER=np.mean([len(z) for z in KEGGS])
    for cluster in KEGGS:
        ALL_KEGGS=ALL_KEGGS+cluster
        
    print record, len(KEGGS), AVERAGE_CLUSTER,len(ALL_KEGGS), len(set(ALL_KEGGS))