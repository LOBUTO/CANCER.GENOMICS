#miRNA project 2
#031514
#Parse genebank chromosomal files

from Bio import SeqIO
from Bio import SeqFeature
import itertools
import subprocess

FILES=subprocess.check_output("ls",cwd="SOFTWARE/pipeline/data/gbs").splitlines()

FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/031514_genomic_coor","w")
FILE_OUT1.write("Chromosome"+"\t"+"Strand"+"\t"+"Start"+"\t"+"End"+"\t"+"miRBase"+"\t"+"gene")

for document in FILES:
    print document
   
    #Need to add +1 in location_start
    #Filter those that are present in mirBase only - those are the only ones that we have targets of
    for record in SeqIO.parse(open(r"SOFTWARE/pipeline/data/gbs/%s"%document), "genbank"):
        FEATURES=filter(lambda x:"ncRNA_class" not in x.qualifiers ,record.features) #Don't want potential mapped processed, only gene coordinates
        FEATURES=filter(lambda w: "db_xref" in w.qualifiers, FEATURES)#Need to have database reference
        FEATURES=filter(lambda y: "miRBase" in list(itertools.chain(*[z.split(":") for z in y.qualifiers["db_xref"]])), FEATURES)
        
        
        for x in FEATURES:
            print x.qualifiers
            mir=[g.split(":")[1] for g in x.qualifiers["db_xref"] if "miRBase" in g][0]
            gene=x.qualifiers['gene'][0]
            
            FILE_OUT1.write("\n"+document.split("chr")[1].split(".")[0]+"\t"+
                            str(x.strand)+"\t"+
                            str(x.location.nofuzzy_start+1)+"\t"+
                            str(x.location.nofuzzy_end)+"\t"+
                            str(mir)+"\t"+
                            str(gene))
                            
FILE_OUT1.close()