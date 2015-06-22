import subprocess
from FUNCTIONS import GENES_TO_UNIPROT_CLASS

GENES=subprocess.check_output("ls", cwd="DATABASES/CANCER_DATA/geneStats").splitlines()
GENES=[X.strip(".stats") for X in GENES]
FILE_OUT=open("DATABASES/CANCER_DATA/CANCER_GENES", "w")
for gene in GENES:
    FILE_OUT.write(gene+"\n")
FILE_OUT.close()