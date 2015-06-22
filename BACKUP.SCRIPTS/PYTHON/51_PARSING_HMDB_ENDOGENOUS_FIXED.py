    #Parsing HMDB_metabolites files to get endogenous metabolites and their protein associations

import subprocess
import xml.etree.ElementTree as ET
def SPLIT_ALL(OPEN_FILE): #Takes file and returns list of appended splitted() content
    DUMMY=[]
    for i in OPEN_FILE:
        DUMMY.append(i.split())
    return DUMMY

#Get list of files

FOLDER_IN=subprocess.check_output("ls", cwd="DATABASES/hmdb_metabolites").splitlines() #40446 metabolites found


#Get endognous metabolites information
FILE_OUT=open("DATABASES/ENDOGENOUS_METABOLITES", "w")
for xml_file in FOLDER_IN:
    FILE=ET.parse("DATABASES/hmdb_metabolites/%s"%xml_file)
    ROOT=FILE.getroot()
    
    #First filter by type of metabolite
    for origin in ROOT.iter("origin"):
        if origin.text=="Endogenous":

            #Get name of metabolite
            NAME=ROOT.find("name").text
            
            #Get ligand PDB id:
            for het_id in ROOT.iter("het_id"):
                HET_ID=str(het_id.text) #Cause if it is empty it produces "Nonetype"
            
            #Get genes associated with ligands
            GENE_NAMES=[]
            for protein_accession in ROOT.iter("protein_accession"):
                GENE_NAMES.append(protein_accession.text)
            
            #Get uniprot ids associated with genes
            UNIPROT=[]
            for uniprot_id in ROOT.iter("uniprot_id"):
                UNIPROT.append(uniprot_id.text)
            
            #Write to file
            FILE_OUT.write(NAME+"|"+HET_ID+" "+"$"+" ")
            for i, j in map(None,GENE_NAMES, UNIPROT):
                FILE_OUT.write(str(i)+"|"+str(j)+" ")
            FILE_OUT.write("\n")

FILE_OUT.close()


#FIX for incorrect het_id codes using SMILES FILES
"""
FILE_IN=open("DATABASES/ENDOGENOUS_METABOLITES")
FILE=FILE_IN.read().splitlines()
FILE_IN.close()

FILE_OUT=open("DATABASES/ENDOGENOUS_METABOLITES_FIXED", "w") 
for i in FILE:
    try:
        print 7
        SMILES_IN=open("Components-smiles-oe.smi") #EITHER SMILES FILE DATABASE FILE HAS IT
        SMILES=SPLIT_ALL(SMILES_IN)
        SMILES_IN.close()        
        if i.split("$")[0].split("|")[1][:-1]!="None" and len(i.split("$")[0].split("|")[1][:-1])>3:
            for j in SMILES:
                if i.split("$")[0].split("|")[0].lower()==j[-1].lower():
                    FILE_OUT.write(j[2] + "|"+j[1]+" "+"$"+" "+i.split("$")[1]+"\n")    
        else:
            FILE_OUT.write(i+"\n")
    
    except IndexError:
        print 8
        SMILES_IN=open("Components-smiles-stereo-cactvs.smi")
        SMILES=SPLIT_ALL(SMILES_IN)
        SMILES_IN.close()    

        if i.split("$")[0].split("|")[1][:-1]!="None" and len(i.split("$")[0].split("|")[1][:-1])>3:
            for j in SMILES:
                if i.split("$")[0].split("|")[0].lower()==j[-1].lower():
                    FILE_OUT.write(j[2] + "|"+j[1]+" "+"$"+" "+i.split("$")[1]+"\n")
        else:
            FILE_OUT.write(i+"\n")
    
FILE_OUT.close()
"""        
        
#FIX for incorrect het_id codes using chem_comp.xml file (name)

#Get NAME-ID list
CHEM_XML_IN=ET.parse("DATABASES/chem_comp.xml")
ROOT=CHEM_XML_IN.getroot()

NAME=ROOT.findall("chemComp/name")
ID=ROOT.findall("chemComp/id")
NAME_ID=[]
for i,j in zip(NAME,ID):
    NAME_ID.append([i.text,j.text])
print len(NAME_ID)
print NAME_ID[0:2]

#FIX
FILE_IN=open("DATABASES/ENDOGENOUS_METABOLITES")
FILE=FILE_IN.read().splitlines()
FILE_IN.close()

FILE_OUT=open("DATABASES/ENDOGENOUS_METABOLITES_FIXED1", "w")
for i in FILE:
    if i.split("$")[0].split("|")[1][:-1]!="None" and len(i.split("$")[0].split("|")[1][:-1])>3:
        for j in NAME_ID:
            if i.split("$")[0].split("|")[0].lower()==j[0].lower():
                FILE_OUT.write(i.split("$")[0].split("|")[0]+"|"+j[1]+" "+"$"+" "+i.split("$")[1]+"\n")
    else:
        FILE_OUT.write(i+"\n")
FILE_OUT.close()