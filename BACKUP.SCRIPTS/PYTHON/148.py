#Cluster run on SPICi with results from Metabolomic network
#04/28/14

import subprocess
"""
#FOR VERSION 1
FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS").splitlines()
FILES=filter(lambda x: x[-1]=="N", FILES)

FOLDER="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS"

for record in FILES:
    FILE_NAME=FOLDER+"/"+record
    subprocess.call("/usr/local/bin/spici -i %s -o %s_5C.cluster -d 0.5 -s 5 -g 0.5 -m 1"%(FILE_NAME,FILE_NAME), shell=True)
"""

"""
#FOR VERSION 2
FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS").splitlines()
FILES=filter(lambda x:"NV2" in x, FILES)

FOLDER="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS"

for record in FILES:
    FILE_NAME=FOLDER+"/"+record
    subprocess.call("/usr/local/bin/spici -i %s -o %s_5C.cluster -d 0.5 -s 5 -g 0.5 -m 1"%(FILE_NAME,FILE_NAME), shell=True)
"""
"""
#FOR VERSION 3
FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS").splitlines()
FILES=filter(lambda x:"NV3" in x, FILES)

FOLDER="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS"

for record in FILES:
    FILE_NAME=FOLDER+"/"+record
    subprocess.call("/usr/local/bin/spici -i %s -o %s_5C.cluster -d 0.5 -s 5 -g 0.5 -m 1"%(FILE_NAME,FILE_NAME), shell=True)

"""

#FOR VERSION X
VERSION="NOMEGA.T.0.5"
FILES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS").splitlines()
FILES=filter(lambda x:VERSION in x, FILES)

FOLDER="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS"

for record in FILES:
    FILE_NAME=FOLDER+"/"+record
    subprocess.call("/usr/local/bin/spici -i %s -o %s_5C.cluster -d 0.5 -s 5 -g 0.5 -m 1"%(FILE_NAME,FILE_NAME), shell=True)