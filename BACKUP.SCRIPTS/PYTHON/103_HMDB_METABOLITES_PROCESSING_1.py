#HMDB METABOLITES PROCESSING 1
#08/14/13

import subprocess, pickle
import xml.etree.ElementTree as ET

#Load files - MODIFY IT BY PLAIN OBJECT LIST PRODUCED IN PROCESSING 0
XML_FILES=subprocess.check_output("ls", cwd="DATABASES/hmdb_metabolites").splitlines()

#Get all metabolites and uniprot/gene as dictionary 
METABOLITE_UNIPROT_DICT={}
for record in XML_FILES:
    print record
    ROOT=ET.parse("DATABASES/hmdb_metabolites/%s"%record).getroot()
    
    #Check if Endogenous to proceed
    ORIGIN=[X.text for X in ROOT.iter("origin")]
    
    if "Endogenous" in ORIGIN and "Drug metabolite" not in ORIGIN:
        #Get name only
        NAME=ROOT.find("name").text.strip()
        
        #Get name synonyms
        SYNONYMS_RECORD=ROOT.find("synonyms")
        NAME=[NAME]+[syn.text for syn in SYNONYMS_RECORD] #ORDER CHANGED
        print NAME
        
        #Get UNIPROTS
        UNIPROT=[X.text.strip() for X in ROOT.iter("uniprot_id")]
        
        #Add to dictionary - Keep only those that have uniprot records
        if len(UNIPROT)>0:
            METABOLITE_UNIPROT_DICT[tuple(NAME)]=UNIPROT #NOOB WARNING - tuple will contain
                                                            #comma after for single elements
        else:
            print "NOPE"

#Store as pickle dictionary
PICKLE_OUT=open("DATABASES/OBJECTS/DICT_METABOLITE_TO_UNIPROTS_HMDB.pi", "w")
pickle.dump(METABOLITE_UNIPROT_DICT, PICKLE_OUT)
PICKLE_OUT.close()
