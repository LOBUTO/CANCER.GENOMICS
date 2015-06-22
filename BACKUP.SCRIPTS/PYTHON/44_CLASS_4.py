# SCRIPT 4 - Make alignments
import subprocess
from Bio import AlignIO


def CLUSTALO_ALIGNIO(FILE_IN, FILE_ALN, FILE_OUT): #ALN FILE and PARSED ALN FILE - Takes in sequences and aligns them using
                                                    #clustal omega
    subprocess.check_output(["CLUSTALO/clustalo", "-i", FILE_IN, "-o", FILE_ALN, "--outfmt=clu", "--force"])
    
    FILE=open(FILE_ALN)
    alignments=AlignIO.parse(FILE, "clustal") #Parses the aligned file to join sequences
    OUTPUT=open(FILE_OUT, "w")
    AlignIO.write(alignments, OUTPUT, "stockholm") #Saves the parsed file in the stockholm format
    
    FILE.close()
    OUTPUT.close()

    #Get rid of "#" and "/" containing lines
    PARSED_IN=open(FILE_OUT)
    PARSED=PARSED_IN.readlines()
    PARSED_IN.close()
    PARSED=filter(lambda x: x[0]!="#" and x[0]!="/", PARSED)
    PARSED_OUT=open(FILE_OUT, "w")
    PARSED_OUT.writelines(PARSED)
    PARSED_OUT.close()
    
SEQ=subprocess.check_output("ls", cwd="CLASS_PROJECT/GROUPS_FILTER_SEQ").splitlines()

for i in SEQ:
    CLUSTALO_ALIGNIO("CLASS_PROJECT/GROUPS_FILTER_SEQ/%s"%i, "CLASS_PROJECT/ALIGNMENT/%s"%i,
                      "CLASS_PROJECT/ALIGNMENT_PARSED/%s"%i)

