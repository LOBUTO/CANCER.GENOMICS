"""
import xml.etree.ElementTree as ET
ROOT=ET.parse("METABOLOMICS/hmdb_metabolites/HMDB00250.xml").getroot()

NAME=ROOT.find("name").text.strip()
UNIPROTS=[X.text.strip() for X in ROOT.findall("protein_associations/protein/uniprot_id")]
GENES=[str(X.text).strip() for X in ROOT.findall("protein_associations/protein/gene_name")]
TYPE=[str(X.text).strip() for X in ROOT.findall("protein_associations/protein/protein_type")]

print NAME
for i in range(len(UNIPROTS)):
    print UNIPROTS[i], GENES[i], TYPE[i]
    
#FOLDER="METABOLOMICS/hmdb_metabolites"
#FILES=os.listdir(FOLDER)[:2]
#print FILES

"""

#Cleaner HMDB2.py
import xml.etree.ElementTree as ET
import os, sys, pickle

HMDB_FOLDER=sys.argv[1]
DRUG_FILE=sys.argv[2]
OUTPUT=sys.argv[3]

#Get xml files from HMDB
FILES=os.listdir(HMDB_FOLDER)
FILES= filter(lambda x: "xml" in x, FILES) #Needed for some mac systems

#Dict output
DICT={}

#STEP 1 - Process files
for record in FILES:
    ROOT=ET.parse(HMDB_FOLDER+"/"+record).getroot()
    
    ORIGIN=[X.text for X in ROOT.iter("origin")]
    ORIGIN_CHECK=set(["Endogenous", "Food"]) | set(ORIGIN)
    
    #Check if origin is from Endogenous or Food source
    if (len(ORIGIN_CHECK)<=2) and (len(ORIGIN)>0):
    
        #Get Name
        NAME=ROOT.find("name").text.strip()        
        
        #Get Synonyms
        SYNONYMS=[SYN.text for SYN in ROOT.findall("synonyms/synonym")]
            
        #Get Smile - KEEP IN MIND THAT SOME MAY HAVE SMILES="None"
        SMILES=str(ROOT.find("smiles").text).strip()
        
        #Get KEGG id - KEEP IN MIND THAT SOME MAY HAVE kegg_id="None"
        KEGG_ID=str(ROOT.find("kegg_id").text).strip()
        
        #Get UNIPROTS
        UNIPROTS=[X.text.strip() for X in ROOT.findall("protein_associations/protein/uniprot_id")]

        #Get Gene names, filter out "None"
        GENES=[str(Y.text).strip() for Y in ROOT.findall("protein_associations/protein/gene_name")]
        GENES=list(set(filter(lambda x: x!="None", GENES)))
 
        #Append to output dict if record has associated genes with it
        if len(GENES)>0:

            DICT[NAME]={"SMILES":"", "GENES":[], "SYNONYMS":[]}
            DICT[NAME]["SMILES"]=SMILES #String
            DICT[NAME]["KEGG_ID"]=KEGG_ID #String
            DICT[NAME]["GENES"]=GENES #List
            DICT[NAME]["SYNONYMS"]=SYNONYMS #List

print len(DICT)

#STEP 2 - Filter for Drugs
FILE_IN1=open(DRUG_FILE)
DRUG_FILTER=[DRUG.upper() for DRUG in FILE_IN1.read().splitlines()]
FILE_IN1.close()

for METABOLITE in DICT.keys():
    if len(set([MET.upper() for MET in [METABOLITE]+DICT[METABOLITE]["SYNONYMS"]]) & set(DRUG_FILTER))>0:
        del DICT[METABOLITE]
print len(DICT)


#STEP 3 - Filter for synonyms
for METABOLITE in DICT.keys():
    print len(DICT)
    
    #Check that we still have the metabolite in the dictionary
    if DICT.has_key(METABOLITE):
        
        #Compare against all other elements in dictionary
        OTHER_METABOLITES=filter(lambda x: x!=METABOLITE, DICT.keys())
        
        for metabolite in OTHER_METABOLITES:    
            
            #Compare all names for similarity
            MAIN_NAMES=set([METABOLITE]+DICT[METABOLITE]["SYNONYMS"])
            COMP_NAMES=set([metabolite]+DICT[metabolite]["SYNONYMS"])
            
            if len(MAIN_NAMES & COMP_NAMES)>0:
                
                #Combine all names into synonym
                DICT[METABOLITE]["SYNONYMS"]=list(set(DICT[METABOLITE]["SYNONYMS"] + list(COMP_NAMES)))
                
                #Combine all genes that bind it
                COMP_GENES=DICT[metabolite]["GENES"]
                DICT[METABOLITE]["GENES"]=list(set(DICT[METABOLITE]["GENES"] + COMP_GENES))
                
                #Delete absorbed entry
                del DICT[metabolite]

print len(DICT)

#STEP 4 - Filter for KEGG_ID Commonoality (ONLY FOR KEGG RELATED ANALYSIS)
# This step may have large differences in SMILES composition which may not be suitable for most analysis
for METABOLITE in DICT.keys():
    print len(DICT)
    
    if DICT.has_key(METABOLITE) and DICT[METABOLITE]["KEGG_ID"]!="None":
        
        #Get its genes
        MAIN_GENES=set(DICT[METABOLITE]["GENES"])
        
        #Make sure meatabolite that are being compared have a KEGG_ID identifier
        OTHER_METABOLITES=filter(lambda y: y!=METABOLITE and DICT[y]["KEGG_ID"]!="None", DICT.keys())
        
        #Compare all genes of metabolites with same KEGG_ID
        #If same KEGG_ID combine only if they have completely overlapping sets of genes
        
        for metabolite in OTHER_METABOLITES:
            
            if DICT[METABOLITE]["KEGG_ID"]==DICT[metabolite]["KEGG_ID"]:
                
                OTHER_GENES=set(DICT[metabolite]["GENES"])
                
                #Only delete if genes are completely overlapping
                if MAIN_GENES==OTHER_GENES:
                    
                    #Combine all names into synonym
                    COMP_NAMES=set([metabolite]+DICT[metabolite]["SYNONYMS"])
                    DICT[METABOLITE]["SYNONYMS"]=list(set(DICT[METABOLITE]["SYNONYMS"] + list(COMP_NAMES)))
                
                    #Combine all genes that bind it
                    DICT[METABOLITE]["GENES"]=list(set(DICT[METABOLITE]["GENES"] + list(OTHER_GENES)))
                    
                    #Delete absorbed entry
                    del (DICT[metabolite])
        
        
#STEP 5 - Remove Water molecule
if ("Water" in DICT)==True:
    del DICT["Water"]
    
#Pickle out
FILE_OUT=open(OUTPUT, "w")
pickle.dump(DICT, FILE_OUT)
FILE_OUT.close()

            
    

