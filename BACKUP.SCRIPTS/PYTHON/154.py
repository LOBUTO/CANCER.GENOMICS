from Bio import SeqIO
import numpy as np

#Load CODON PROBABILITY table
FILE_IN1=open("DATABASES/CODON.PROBABILITY")
COD_PROB=[x.split("\t") for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

#Build SILENT and NON_SILENT codon probability dictionaries
SILENT_DICT=dict((x[0], x[2]) for x in COD_PROB)
NON_SILENT_DICT=dict((x[0], x[1]) for x in COD_PROB)

FILE_OUT=open("DATABASES/UNIPROT/073114_HUMAN_AA_SILENT_RATIO","w")
FILE_OUT.write("Hugo_Symbol"+"\t"+"SILENT.COUNT"+"\t"+"NON.SILENT.COUNT"+"\t"+"SNS_SUM_RATIO"+
               "\t"+"SILENT_MEAN"+"\t"+"NON_SILENT_MEAN"+"\t"+"SNS_MEAN_RATIO"+
               "\t"+"SILENT_MEDIAN"+"\t"+"NON_SILENT_MEDIAN"+"\t"+"SNS_MEDIAN_RATIO"+
               "\t"+"SEQUENCE_LENGTH")

for record in SeqIO.parse("DATABASES/UNIPROT/073014_HUMAN_AA_REVIEWED.fasta", "fasta"):
    if "GN=" in record.description:
        print record.id, record.description.split("GN=")[1].split()[0], record.seq
        GENE=str(record.description.split("GN=")[1].split()[0])
        SEQUENCE=filter(lambda y: y in SILENT_DICT , [z for z in str(record.seq)]) #Make sure that we have amino acids in our dicts
        SEQUENCE_LENGTH=len(SEQUENCE)
        
        SILENT_VECTOR=[float(SILENT_DICT[x]) for x in SEQUENCE]
        NON_SILENT_VECTOR=[float(NON_SILENT_DICT[x]) for x in SEQUENCE]
        
        SILENT_SUM=sum(SILENT_VECTOR)
        NON_SILENT_SUM=sum(NON_SILENT_VECTOR)
        SILENT_MEAN=np.mean(SILENT_VECTOR)
        NON_SILENT_MEAN=np.mean(NON_SILENT_VECTOR)
        SILENT_MEDIAN=np.median(SILENT_VECTOR)
        NON_SILENT_MEDIAN=np.median(NON_SILENT_VECTOR)
        
        FILE_OUT.write("\n"+GENE+"\t"+str(SILENT_SUM)+"\t"+str(NON_SILENT_SUM)+"\t"+str(SILENT_SUM/NON_SILENT_SUM)+
                       "\t"+str(SILENT_MEAN)+"\t"+str(NON_SILENT_MEAN)+ "\t"+str(SILENT_MEAN/NON_SILENT_MEAN)+
                       "\t"+str(SILENT_MEDIAN)+"\t"+str(NON_SILENT_MEDIAN)+ "\t"+ str(SILENT_MEDIAN/NON_SILENT_MEDIAN)+
                       "\t"+str(SEQUENCE_LENGTH))

FILE_OUT.close()