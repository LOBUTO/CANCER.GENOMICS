#miRNA project 1
#031314
#Parse miRNA EMBL

from Bio import SeqIO
from Bio import SeqFeature
import itertools

FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/031514_ACC_SEQ","w")
FILE_OUT1.write("miRBase"+"\t"+"Sequence"+"\t"+"Mature"+"\t"+"Start"+"\t"+"End")

for record in SeqIO.parse(open(r"DATABASES/CANCER_DATA/miRNA/miRNA.dat.dat"), "embl"):
    ID=record.id
    MATURE_LOCATION=[[X.location.nofuzzy_start+1, X.location.nofuzzy_end] for X in record.features]
    PRODUCT= [X[1][0] for Y in record.features for X in Y.qualifiers.items() if X[0]=="product"]
    
    MATURE_LOCATION=filter(lambda y:y[0]!=y[1], MATURE_LOCATION)
    
    if "hsa" in list(itertools.chain(*[z.split("-") for z in PRODUCT])): #Filter for human (hsa)
        
        for index in list(xrange(len(MATURE_LOCATION))):
            
            FILE_OUT1.write("\n"+str(ID)+"\t"+str(record.seq)+"\t"+str(PRODUCT[index])+"\t"+
                            str(MATURE_LOCATION[index][0])+"\t"+str(MATURE_LOCATION[index][1]))

FILE_OUT1.close()
