#EMBOSS_ALIGNIO FUNCTION
#8/4/13
#Dependencies:
#    Biopython
#Takes in two FASTA!!! files, performs global pairwise alignment and parses the files to FILE_OUT ready for asterisk replacement method

def EMBOSS_ALIGNIO(FASTA1, FASTA2, FILE_OUT, LOGS):    
    import subprocess
    from Bio import AlignIO
    
    #Do global alignment
    subprocess.check_output(["/usr/local/bin/needle", FASTA1, FASTA2, 
                      "gapopen=11", "gapextend=1", "%s/EA1"%LOGS])
    
    #Convert to stockholm format for easier parsing
    FILE_IN=open("%s/EA1"%LOGS)
    ALIGNIO=AlignIO.parse(FILE_IN, "emboss")
    OUTPUT=open("%s/EA2"%LOGS, "w")
    AlignIO.write(ALIGNIO, OUTPUT, "stockholm")
    
    FILE_IN.close()
    OUTPUT.close()
    
    #Get rid of "#" and "/" containing lines
    PARSED_IN=open("%s/EA2"%LOGS)
    PARSED=PARSED_IN.readlines()
    PARSED_IN.close()
    PARSED=filter(lambda x: x[0]!="#" and x[0]!="/", PARSED)
    PARSED_OUT=open(FILE_OUT, "w")
    PARSED_OUT.writelines(PARSED)
    PARSED_OUT.close()    

#TEST
#EMBOSS_ALIGNIO("LOGS/LOG1.fasta", "LOGS/LOG2.fasta", "LOGS/OUT.PARSED", "LOGS")
EMBOSS_ALIGNIO("LOGS/PRTUR1.fasta", "LOGS/PRTUR2.fasta", "LOGS/TESEE", "LOGS")
