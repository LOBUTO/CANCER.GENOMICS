#LEARNING EMBOSS NEEDLE, EMBOSS WATER, BLAST STANDALONE AND BIOPYTHON BLAST
import subprocess
from Bio import AlignIO

#Water for local alignment and needle for global alignment
"""
FILE_OUT=subprocess.Popen(["/usr/local/bin/water", "LOGS/LOG1.fasta", "LOGS/LOG2.fasta", 
                  "gapopen=11", "gapextend=1", 
                  "-stdout", "stdout"], stdout=subprocess.PIPE)
ALN=FILE_OUT.communicate()[0]
print ALN
"""
#print 7
#NEW=AlignIO.parse(FILE_OUT.stdout, "emboss")
#for row in NEW: print row


#BLAST FROM STANDALONE
#Running a query (QUERY IS UNIPROT AND SUBJECT IS PDB)
"""
import subprocess
#To get complete alignment
RUN1=subprocess.Popen(["-c", "blastp -query=LOGS/LOG1.fasta -subject=LOGS/LOG2.fasta"], stdout=subprocess.PIPE, shell=True)
PARSED1=RUN1.communicate()[0]
print PARSED1

#To get specific data out of it use format 7
RUN2=subprocess.Popen(["-c", "blastp -query=LOGS/LOG1.fasta -subject=LOGS/LOG2.fasta -outfmt '7  pident qcovs length'"],
                      stdout=subprocess.PIPE, shell=True)
PARSED2=RUN2.communicate()[0]
print PARSED2
"""
#BLAST FROM BIOPYTHON - SLOWER THAN STANDALONE
from Bio.Blast import NCBIWWW, NCBIXML
FASTA_STRING=open("LOGS/LOG4.fasta").read()
RESULTS=NCBIWWW.qblast("blastp", "swissprot", FASTA_STRING)
print RESULTS
RECORD=NCBIXML.read(RESULTS)
print RECORD

for alignment in RECORD.alignments:
    print "START"
    print alignment
    print 7
    for hsp in alignment.hsps:
        print hsp
