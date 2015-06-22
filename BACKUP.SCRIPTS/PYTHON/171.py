###Parse RECON xml file####
import re

def anno_class (anno_list):
    CHEBI=[anno.split("CHEBI:")[1] for anno in anno_list if "chebi" in anno] 
    KEGG=[anno.split("/")[-1] for anno in anno_list if "kegg" in anno]
    HMDB=[anno.split("hmdb/")[1] for anno in anno_list if "hmdb" in anno]
    return {"CHEBI":CHEBI, "KEGG":KEGG, "HMDB":HMDB}

def body_class (body_list):
    SUBSYSTEM=filter(lambda x: "SUBSYSTEM" in x, body_list)[0].split("SUBSYSTEM:")[1].strip()
    if SUBSYSTEM=="":
        SUBSYSTEM="NONE"
    EC=filter(lambda x: "EC" in x, body_list)[0].split("Number:")[1].strip()
    if EC=="":
        EC="NONE"
    return {"SUBSYSTEM":SUBSYSTEM, "EC":EC}

def obo_parser (obo_file):
    
    #Start dictionary and mock entry
    DICT={}
    ID="MOCK"
    KEGG_ID=""
    ALT=[]
    TYPE_DEF_COUNT=0
    
    with open(obo_file, "r") as infile:
        for line in infile:
            line=line.strip()
            
            #Process only terms
            if line=="[Term]":
                
                #Only add to dict if we haven't encountered "type defs"
                if TYPE_DEF_COUNT==0 and KEGG_ID!="":
                    
                    #Add to dict only if we have KEGG_IDs and one per each main name and synonym
                    for chebi in [ID]+ALT:
                        DICT[chebi]=KEGG_ID
                    
                    ID=""
                    KEGG_ID=""
                    ALT=[]
                else:
                    ID=""
                    KEGG_ID=""
                    ALT=[]
                    
                    #And reset type_def_count
                    TYPE_DEF_COUNT=0
                
            #Obtain main id
            elif line[0:10]=="id: CHEBI:":
                ID=line.split("CHEBI:")[-1].strip()
            
            #Obtain alternative if any 
            elif line[0:14]=="alt_id: CHEBI:":
                ALT=ALT+[line.split("CHEBI:")[-1].strip()]
            
            #Obtain KEGG_ID if any
            elif line[0:12]=="xref: ChEBI:" and "KEGG COMPOUND" in line and line.split()[1].split(":")[1][0:2]=="C0":
                KEGG_ID=line.split()[1].split(":")[1]
            
            elif line[0:20]=="xref: KEGG COMPOUND:" and line.split()[2].split(":")[1][0:2]=="C0":
                KEGG_ID=line.split()[2].split(":")[1]
            
            #Watch out for "type_defs"
            elif line=="[Typedef]":
                TYPE_DEF_COUNT=TYPE_DEF_COUNT=0
    
    #Delete mock entry
    if "MOCK" in DICT:
        DICT.pop("MOCK")
    return(DICT)

def weight_dict (periodic_file):
    
    file_in=open(periodic_file)
    ch_dict=dict((x.split("=")[0].strip(), x.split("=")[1].strip()) for x in file_in.read().splitlines() )
    file_in.close()
    
    return ch_dict

print "EXECUTING"

##Obtain weight_dict
el_dict=weight_dict("/Users/jzamalloa/Desktop/FOLDER/GIT/periodic/periodic/mass.py")

##Obtain CHEBI dict
OBO_DICT=obo_parser("DATABASES/CHEBI/chebi.obo.txt")

##Open files to write to
FILE_OUT1=open("DATABASES/RECON/042215.PROCESSED.METABOLITES", "w")
FILE_OUT1.write("RECON.ID"+ "\t"+ "NAME"+ "\t"+ "COMPARTMENT"+"\t"+ "KEGG_ID"+"\t"+ "HMDB"+"\t"+ "CHEBI"+"\t"+"EHMN"+ "\t"+"WEIGHT")
FILE_OUT2=open("DATABASES/RECON/042215.PROCESSED.GENE", "w")
FILE_OUT2.write("MODIFIERS"+"\t"+ "Hugo_Symbol"+"\t"+"COMPARTMENT")

#import xml.etree.ElementTree as ET
from lxml import etree

##To parse the file so you could start iterating and using it as a tree (nested lists)
FILE=etree.parse("DATABASES/RECON/recon2_02/recon2.v02.xml")

ROOT=FILE.getroot()

#Get the xlmns
xlmns=ROOT.tag.rstrip("sbml")

####First annotate SPECIES####
for child in ROOT[0][4]:
    
    #Make sure we have a name to work with 
    if "name" in child.attrib.keys():
        
        #Obtain identifiers
        ID=child.attrib["id"]
        NAME=child.attrib["name"]
        COMPARTMENT=child.attrib["compartment"]
    
        #Check if we are dealing with metabolite...
        if "charge" in child.attrib.keys():
    
            #print child.attrib["id"], child.attrib["name"], child.attrib["compartment"] #take care of "name":"null"
            #Find external annotations, ie. KEGG, HDBM, CHEBI
            ANNO=[dict(anno.attrib).values()[0] for anno in child.iter("*") if len(anno.attrib)>0]
            
            #Only annotation if annotation exist
            if len (ANNO)>0:
                ANNO=anno_class(ANNO)
            
            else: 
                #Keep in mind that this is for metabolites only
                ANNO={"CHEBI":[], "KEGG":[], "HMDB":[]}
            
            #Add EHMN annotation if present
            ANNO["EHMN"]=[]
            if len([x.tag for x in child if "notes" in x.tag])>0:
                EHMN=[y.text.split(":")[1].strip() for y in child[0][0] if "EHMN" in y.text] 
                ANNO["EHMN"]=EHMN
            
            #Use EHMN for KEGG annotation if posseses the correct ID format
            if ANNO["KEGG"]==[] and ANNO["EHMN"]!=[] and ANNO["EHMN"][0][0:2]=="C0":
                ANNO["KEGG"]=ANNO["EHMN"] 
            
            #Use CHEBI dict if KEGG still not present
            elif ANNO["KEGG"]==[] and ANNO["CHEBI"]!=[] and ANNO["CHEBI"][0] in OBO_DICT:
                ANNO["KEGG"]=[OBO_DICT[ANNO["CHEBI"][0]]]
            
            ###MANUALLY CORRECT ANNOTATION#####
            if ID in ["M_cl_e", "M_cl_c", "M_cl_l", "M_cl_n", "M_cl_m", "M_cl_g", "M_cl_x", "M_cl_r"]:
                ANNO["KEGG"]=["C00698"]
                
            ###################################
            
            #Finally annotate formula
            FORMULA=[x.text for x in child[0][0] if "FORMULA" in x.text][0].split(":")[1].strip()
            if len(FORMULA)>0 and "R" not in FORMULA  and "X" not in FORMULA:
                
                FORMULA=FORMULA.rstrip("X")
                
                #Split formula into element and integer bases
                FORMULA=re.findall('[A-Z][^A-Z]*', FORMULA)
                print (ID,FORMULA)
                CUM_SUM=0
                for sub in FORMULA:
                    EL="".join([x for x in sub if x.isdigit()==False])
                    DIGIT="".join([x for x in sub if x.isdigit()])
                    if DIGIT=="":
                        DIGIT=1
                    else:
                        DIGIT=int(DIGIT)
                    CUM_SUM=CUM_SUM + float(el_dict[EL])*DIGIT 
                    
            else:
                CUM_SUM="NONE"
            
            #Write to file
            FILE_OUT1.write("\n"+ID+"\t"+NAME+"\t"+COMPARTMENT)
            for anno in ["KEGG","HMDB","CHEBI", "EHMN"]:
                if len(ANNO[anno])!=0:
                    FILE_OUT1.write("\t"+ANNO[anno][0])
                else:
                    FILE_OUT1.write("\t"+"NONE")
            FILE_OUT1.write("\t"+str(CUM_SUM))
            
        #OR gene....
        else:
            #Write to file
            for name in NAME.split(":"):
                FILE_OUT2.write("\n"+ID+"\t"+name+"\t"+COMPARTMENT)
        
#Close files
FILE_OUT1.close()
FILE_OUT2.close()

####Second pass at annotating metabolites####
FILE_OUT1=open("DATABASES/RECON/042215.PROCESSED.METABOLITES")
MET_TABLE=[x.split("\t") for x in FILE_OUT1.read().splitlines()[1:]]
MET_DICT=dict((x[0][0:-1],x[3]) for x in MET_TABLE if x[3]!="NONE")
FILE_OUT1.close()
print MET_TABLE[1:5]
print MET_DICT.keys()[1:5]
print MET_DICT["M_glutcoa_"]

FILE_OUT1=open("DATABASES/RECON/042215.PROCESSED.METABOLITES", "w")
FILE_OUT1.write("RECON.ID"+ "\t"+ "NAME"+ "\t"+ "COMPARTMENT"+ "\t"+ "KEGG_ID"+"\t"+ "HMDB"+"\t"+ "CHEBI"+"\t"+"EHMN"+ "\t"+"WEIGHT")

for line in MET_TABLE:
    
    if line[3]=="NONE" and line[0][:-1] in MET_DICT:
        FILE_OUT1.write("\n"+ "\t".join(line[0:3]) + "\t"+MET_DICT[line[0][:-1]] + "\t" + "\t".join(line[4:]) ) 
    else:
        FILE_OUT1.write("\n" + "\t".join(line))
    
FILE_OUT1.close()

####Then annotate REACTIONS#####
#Open file to write reactions
FILE_OUT4=open("DATABASES/RECON/042215.PROCESSED.REACTIONS", "w")
FILE_OUT4.write("REACTION.ID"+ "\t"+ "REACTION.NAME"+"\t"+"REVERSIBLE"+"\t"+"SUBSYSTEM"+"\t"+"EC"+
                "\t"+"SUBSTRATES"+"\t"+"PRODUCTS"+"\t"+"MODIFIERS")

print ("REACTIONS")
for child in ROOT[0][5]:
    
    #Obtain identifiers
    ID=child.attrib["id"]
    NAME=child.attrib["name"]
    REVERSIBLE=child.attrib["reversible"]
    
    #Get subsystem (proxy of pathway) and EC
    BODY=[body.text for body in child.findall("".join([xlmns,"notes"]))[0][0]]
    BODY=body_class(BODY)
    
    #Get reactants and products if they exist
    if len(child.findall("".join([xlmns, "listOfReactants"])))>0:
        REACTANTS=[s.attrib["species"] for s in  child.findall("".join([xlmns, "listOfReactants"]))[0]]
        
    else:
        REACTANTS=["NONE"]
    
    if len(child.findall("".join([xlmns, "listOfProducts"])))>0:
        PRODUCTS=[s.attrib["species"] for s in  child.findall("".join([xlmns, "listOfProducts"]))[0]]
    else:
        PRODUCTS=["NONE"]
    
    #Get modifier genes if they exist
    if len(child.findall("".join([xlmns, "listOfModifiers"])))>0:
        MODIFIERS=[dict(s.attrib).values()[0] for s in  child.findall("".join([xlmns, "listOfModifiers"]))[0]]
    else:
        MODIFIERS=[]
    
    #Write to file if gene associated with it
    if len(MODIFIERS)>0:
        
        for mod in MODIFIERS:
            FILE_OUT4.write("\n"+ID+"\t"+NAME+"\t"+REVERSIBLE+"\t"+BODY["SUBSYSTEM"]+"\t"+BODY["EC"]+"\t"+
                        "|".join(REACTANTS)+"\t"+"|".join(PRODUCTS)+"\t"+mod)
            
#Close file
FILE_OUT4.close()
