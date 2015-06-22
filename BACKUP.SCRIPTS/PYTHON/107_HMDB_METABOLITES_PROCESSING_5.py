#HMDB METABOLITES PROCESSING 5 (USED FOR 4.1) - HELPER
#08/15/13
#Create dictionary object of metabolites to smiles codes
#USED IN THE CONTRUCTION OF STEP 4.1

import subprocess, pickle
import xml.etree.ElementTree as ET

#Get list of xml files
XML_FILES=subprocess.check_output("ls", cwd="DATABASES/hmdb_metabolites").splitlines()

#Get names and metabolites
FILE_OUT1=open("DATABASES/HMDB_SOURCE/HAVE_SMILES.smi","w")
FILE_OUT2=open("DATABASES/HMDB_SOURCE/HAVE_NO_SMILES.smi","w")
MET_TO_SMILES_DICT={}
for record in XML_FILES:
    print record
    ROOT=ET.parse("DATABASES/hmdb_metabolites/%s"%record).getroot()
    
    #Get smiles first
    SMILES=ROOT.find("smiles").text
    
    #Get name second
    NAME=ROOT.find("name").text.strip()
    
    #Store in file
    try:
        if str(SMILES)!="None":
            FILE_OUT1.write(SMILES+"\t"+NAME+"\n")
            MET_TO_SMILES_DICT[NAME]=SMILES
        else:
            FILE_OUT2.write(NAME+"\n")
    
    except UnicodeEncodeError:
        NAME=NAME.encode(encoding="UTF-8") #Because some have unicode encoding like +/- in one line        
        if str(SMILES)!="None":
            FILE_OUT1.write(SMILES+"\t"+NAME+"\n")
            MET_TO_SMILES_DICT[NAME]=SMILES
        else:
            FILE_OUT2.write(NAME+"\n")       
        
FILE_OUT1.close()
FILE_OUT2.close()
PICKLE_OUT1=open("DATABASES/HMDB_SOURCE/OBJECTS/DICT_METABOLITE_TO_SMILES.pi", "w")
pickle.dump(MET_TO_SMILES_DICT, PICKLE_OUT1)
PICKLE_OUT1.close()
