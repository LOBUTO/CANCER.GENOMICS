#SCRIPT TO CALL DARIO'S TCGA PIPELINE
#12/17/13
#Remeber to modify pipeline so that refseq folder is manually called, otherwise this script won't work

import subprocess

#Get all cancer names
CANCERS=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS").splitlines()

#Get .maf files
MAF_FILES=[[Y,X] for Y in CANCERS for X in subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/"+Y).splitlines() if "maf" in X]

#Run Dario's script
for record in MAF_FILES:
    #.maf file location
    MAF_FILE="DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/"+record[0]+"/"+record[1]
    
    #Create cancer directory in mutDir
    subprocess.call("mkdir SOFTWARE/pipeline/step.01/mutDir/%s"%record[0], shell=True)
    
    #Call function
    subprocess.call("python SOFTWARE/pipeline/code/python/processTCGAMut.py %s SOFTWARE/pipeline/data/hsRef SOFTWARE/pipeline/data/gbs \
    SOFTWARE/pipeline/data/human.protein.faa SOFTWARE/pipeline/step.01/refSeqs SOFTWARE/pipeline/step.01/mutDir/%s \
    DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/%s.mut"%(MAF_FILE, record[0], record[0]), shell=True)
    
    
