######HMDB_PRE_PROCESS.py##########
#061114
#- Process pickle output of HMDB2.py to further identify KEGG ids with output from
#  KEGG_ENZYMES.py.
#- Converts pickle to table and removes any metabolite that does not contain any KEGG ID
#Based on 3.py(PC), 032914

import pickle, sys

HMDB_PICKLE=sys.argv[1]
KEGG_FILE=sys.argv[2]
OUTPUT_FILE=sys.argv[3]

#Make KEGG_ID: [Metabolites] dictionary
FILE_IN1=open(KEGG_FILE) #Straight from KEGG_ENZYMES.py
KEGG=[x.split("\t") for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

KEGG_ID_DICT=dict((x[1], []) for x in KEGG) #KEGG_ID empty dictionary

for record in KEGG: #KEGG_ID: [Products]
    KEGG_ID_DICT[record[1]]=list(set(KEGG_ID_DICT[record[1]]+[record[2]]))

#Apply to those HMDB metabolites which don't have yet a KEGG_ID - Give them KEGG_ID if not have it already
FILE_IN2=open(HMDB_PICKLE)#STRAIGHT FROM HMDB2 function
DICT=pickle.load(FILE_IN2)
FILE_IN2.close()
print "TOTAL",len(DICT)
print "WITH KEGG",len(filter(lambda x: DICT[x]["KEGG_ID"]!="None", DICT.keys())) #Initially those that have a KEGG_ID

for METABOLITE in DICT.keys():
    
    #For those that don't have an assigned KEGG_ID
    if DICT[METABOLITE]["KEGG_ID"]=="None":
        
        METABOLITE_NAMES=[METABOLITE]+DICT[METABOLITE]["SYNONYMS"]
        
        for KEGG_ID in KEGG_ID_DICT.keys():
            if len(set(METABOLITE_NAMES) & set(KEGG_ID_DICT[KEGG_ID]))>0:
                
                #If common names are found then assing KEGG_ID
                DICT[METABOLITE]["KEGG_ID"]=KEGG_ID

print "TOTAL", len(DICT)                
print "WITH KEGG",len(filter(lambda x: DICT[x]["KEGG_ID"]!="None", DICT.keys())) #Initially those that have a KEGG_ID

#Finally, re-filter for KEGG commonality (Consolidating within dictionary)
#Merging entries by overlapping genes
for METABOLITE in DICT.keys():
    
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

print "TOTAL", len(DICT)                
print "WITH KEGG",len(filter(lambda x: DICT[x]["KEGG_ID"]!="None", DICT.keys()))

#Write TABLE 2 - BASED ON FIRST NAME, but include KEGG identifier, for reference since there are duplicates
#Keep in mind that this table only has those that have KEGG identifiers!!
FILE_OUT1=open(OUTPUT_FILE, "w")
FILE_OUT1.write("METABOLITE"+"\t"+"KEGG_ID"+"\t"+"GENE")
for record in DICT.keys():
    
    #Only write to Table 2 if it has a KEGG_ID!!!!!
    if DICT[record]["KEGG_ID"]!="None":
        for gene in list(set(DICT[record]["GENES"])): #Make sure to unique gene sets just in case
            FILE_OUT1.write("\n"+record+"\t"+DICT[record]["KEGG_ID"]+"\t"+gene)
        
FILE_OUT1.close()

