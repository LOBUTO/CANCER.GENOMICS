#Prepare file to be parsed for networkx

FILE=open("NETWORK/METABOLITES_ENDOGENOUS_FIXED").read().splitlines()

FILE_OUT=open("NETWORK/ENDOGENOUS.ADJ", "w")
for line in FILE:
    #Get ligand ID
    LIG=line.split("$")[0].split("|")[0]
    LIG="_".join(LIG.split()) #IMPORTANT - JOINING SEPARATE WORDS IN LIGAND BY "_"
    
    #Get gene IDs
    GENES=[]
    for gene in line.split("$")[1].split():
        GENES.append(gene)
    
    #Write it per line
    FILE_OUT.write(u"%s"%LIG.encode("utf-8")+" ") #NETWORKX READS UTF-8 FORMAT SO NEED TO CONVERT IT
    
    for ID in GENES:
        FILE_OUT.write(u"%s"%ID.encode("utf-8")+" ")
    
    FILE_OUT.write("\n")

FILE_OUT.close()

    
        