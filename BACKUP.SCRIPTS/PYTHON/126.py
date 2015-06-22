#01/18/14

def LOC_REPLACE(ITEMS_LIST, REPLACE_DICT):
    
    #Check if it is in list first
    IN_LIST=filter(lambda x:  x in ITEMS_LIST , REPLACE_DICT.keys())
    
    if len(IN_LIST)>0:
        for item in IN_LIST:
            ITEMS_LIST[ITEMS_LIST.index(item)]=REPLACE_DICT[item]
        return ITEMS_LIST
    else:
        return ITEMS_LIST

FILE_IN2=open("DATABASES/CANCER_DATA/THPA/subcellular_location.csv")
LOC=[ [Y.strip('"') for Y in X.split(",")] for X in FILE_IN2.read().splitlines()[1:]]
FILE_IN2.close()

##NOTES:
#Replace:
#    'Nucleus but not nucleoli' with "Nucleus"
#    "Cytoskeleton (____)" with "Cytoskeleton"
#    "Nucleoli" with "Nucleolus"
#    "Plasma membrane" with "Cell membrane"
#    "Mitochondria" with "Mitochondrion"
#    "Focal Adhesions" with "Focal adhesion"
#    "Cell Junctions" with "Cell junction"
#    "Nuclear membrane" with "Nucleus membrane"
#For the time being "Main location" and "Other location" will be combined

FILE_OUT1=open("DATABASES/CANCER_DATA/THPA/THPA_LOCATION","w")
FILE_OUT1.write("GENE"+"\t"+"LOCATION")

for record in LOC:
    
    ENSEMBLE=record[0]
    
    if len(record[1])>0:
        MAIN=record[1].split(";")
    else:
        MAIN=[]
    
    if len(record[2])>0:
        OTHER=record[2].split(";")
    else:
        OTHER=[]
    
    NAMES=list(set(MAIN+OTHER))
    if len(NAMES)>0:
        
        REPLACE_DICT={"Cytoskeleton (Microtubules)":"Cytoskeleton", "Nucleus but not nucleoli":"Nucleus",
                      "Cytoskeleton (Cytokinetic bridge)":"Cytoskeleton", "Cytoskeleton (Actin filaments)":"Cytoskeleton",
                      "Cytoskeleton (Intermediate filaments)":"Cytoskeleton", "Nucleoli":"Nucleolus", 
                      "Plasma membrane":"Cell membrane", "Mitochondria":"Mitochondrion", "Focal Adhesions":"Focal adhesion",
                      "Cell Junctions":"Cell junction", "Nuclear membrane": "Nucleus membrane"}
        NAMES=LOC_REPLACE(NAMES, REPLACE_DICT)
        
        for location in NAMES:
            FILE_OUT1.write("\n"+ENSEMBLE+"\t"+location)

FILE_OUT1.close()
        