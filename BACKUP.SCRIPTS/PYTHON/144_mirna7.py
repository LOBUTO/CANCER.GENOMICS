#miRNA project 7
#031814
#Prepare miRNA hairpin fasta file by extracting only human sequences

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#Filter for human (hsa)
RECORD=[record for record in SeqIO.parse(open(r"DATABASES/CANCER_DATA/miRNA/hairpin.fa.fa"), "fasta") if "hsa" in record.id]

NEW_RECORD=[]
for record in RECORD:
    SEQ=str(record.seq.back_transcribe())
    
    NEW_REC=SeqRecord(Seq(SEQ), id=record.id, name=record.name, description=record.description)
    NEW_RECORD.append(NEW_REC)
    

#Write out filtered fasta
FILE_OUT1=open("DATABASES/CANCER_DATA/miRNA/hairpin.hsa.fa", "w")
SeqIO.write(NEW_RECORD, FILE_OUT1, "fasta")
FILE_OUT1.close()

#mirna fams give better coverage