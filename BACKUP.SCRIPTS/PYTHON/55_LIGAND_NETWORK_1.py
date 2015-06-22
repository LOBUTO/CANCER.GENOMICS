#PARSING HMDB METABOLITES - FIX FOR SYNONYMS
import subprocess
import xml.etree.ElementTree as ET

FILE_IN=subprocess.check_output("ls", cwd="DATABASES/hmdb_metabolites").splitlines()
print len(FILE_IN)

#First get metabolite names that are not synonyms
MB_UNIQ=open("NETWORK/MB_UNIQ", "w")

SYN_LIST=[]
for item in FILE_IN:
    print item
    FILE=ET.parse("DATABASES/hmdb_metabolites/%s"%item)
    ROOT=FILE.getroot()
    
    #First check that it is Endogenous
    for origin in ROOT.iter("origin"):
        if origin.text=="Endogenous":
            
            #Get synonyms of item
            SYN=[]
            for SYNO in ROOT.iter("synonym"):
                SYN.append(SYNO.text)
            
            #Check if has any synonyms at all first
            if len(SYN)>0:
            
                #Get prot_acc
                PROT_LIST=[]
                for prot in ROOT.iter("protein_accession"):
                    PROT_LIST.append(prot.text)
                PROT_LIST=sorted(PROT_LIST)
                
                PL=""
                for s in PROT_LIST:
                    PL=PL+s+"_"
                PROT_LIST=PL
                                
                #First check if synonyms exist
                CROSS_SYN=[]
                X=0
                for i in SYN:
                    if any(i==j.split("|")[0] for j in SYN_LIST):
                        X=X+1 #If synonym does exist add 1 to X
                        CROSS_SYN.append(i)
                    
                #If term has no synonym elsewhere then write
                if X==0:
                    MB_UNIQ.write(item+" "+ROOT.find("name").text+"\n")
                
                #However, if it does exist check if all other terms that contain that synonym have different protein sets
                if X>0:
                    Y=0
                    for l in CROSS_SYN:
                        for ll in SYN_LIST:
                            if l==ll.split("|")[0] and ll.split("|")[2]==PROT_LIST: #If the synonym in our file is equal to a synonym in another file
                                                                                    #and has the same protein set then add Y+1
                                Y=Y+1
                    
                    if Y==0: #If the synonym does not have another synonym in another file with the same protein set then write it
                        MB_UNIQ.write(item+" "+ROOT.find("name").text+"\n")
                
                for p in SYN:
                    SYN_LIST.append(p+"|"+item+"|"+PROT_LIST)
            
            #If it doesn't contain any synonyms then just write it:
            else:
                MB_UNIQ.write(item+" "+ROOT.find("name").text+"\n")
                    
MB_UNIQ.close()

