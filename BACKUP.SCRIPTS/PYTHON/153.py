#NETWORK SPICI ANALYSIS 
#071514
#Replace 148.py

import subprocess

#Load *.Table.W.V
FILE_IN1=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/071514.BRCA.Table.W.V")
TABLES=[x.split() for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

#Get BETAS
BETAS=list(set([x[3] for x in TABLES]))

#Construct table per beta and get clusters
for beta in BETAS:
    print beta
    #Get individual beta.table
    TABLES_BETA=filter(lambda x: x[3]==beta, TABLES)
    
    #Write to table - DON'T WRITE 0 EDGES!!
    FILE_OUT=open("PIPELINES/METABOLIC.DRIVERS/CLUSTER.ANALYSIS/BRCA/NETWORKS/071514_"+beta,"w")
    for record in TABLES_BETA:
        if record[2]!="0":
            FILE_OUT.write(record[0]+"\t"+record[1]+"\t"+record[2]+"\n")
    FILE_OUT.close()
    
#Get cluster for each file
FILES=subprocess.check_output("ls", cwd="PIPELINES/METABOLIC.DRIVERS/CLUSTER.ANALYSIS/BRCA/NETWORKS").splitlines()

FOLDER_IN="PIPELINES/METABOLIC.DRIVERS/CLUSTER.ANALYSIS/BRCA/NETWORKS"
FOLDER_OUT="PIPELINES/METABOLIC.DRIVERS/CLUSTER.ANALYSIS/BRCA/SPICI.RESULTS"

for record in FILES:
    FILE_NAME_IN=FOLDER_IN+"/"+record
    FILE_NAME_OUT=FOLDER_OUT+"/"+record
    subprocess.call("/usr/local/bin/spici -i %s -o %s.cluster -d 0.7 -s 3 -g 0.5 -m 0"%(FILE_NAME_IN,FILE_NAME_OUT), shell=True)

