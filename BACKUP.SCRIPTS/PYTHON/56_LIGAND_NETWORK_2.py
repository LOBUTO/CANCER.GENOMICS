#BUILD NON-SYNONYMS FILE FOR NETWORKX
import xml.etree.ElementTree as ET

UNIQ_XML=open("NETWORK/MB_UNIQ").read().splitlines()

FILE_OUT=open("NETWORK/METABOLITES_ENDOGENOUS", "w")
for item in UNIQ_XML:
    XML=item.split()[0]    
    ROOT=ET.parse("DATABASES/hmdb_metabolites/%s"%XML).getroot()
    
    #Get name
    NAME=ROOT.find("name").text
    
    #Get het_id
    for het_id in ROOT.iter("het_id"):
        HET_ID=str(het_id.text)
        
    #Get HMDB protein_accession IDs associated with ligand
    PROT_ACC=[]
    for protein_accession in ROOT.iter("protein_accession"):
        PROT_ACC.append(protein_accession.text)
    
    #Get Uniprot associated with it
    UNIPROT=[]
    for uniprot_id in ROOT.iter("uniprot_id"):
        UNIPROT.append(uniprot_id.text)
    
    if len(PROT_ACC)!=len(UNIPROT):print "ERROR"
    
    #Write to file
    FILE_OUT.write(NAME+"|"+HET_ID+"$")
    for i, j in map(None,PROT_ACC, UNIPROT):
        FILE_OUT.write(str(i)+"|"+str(j)+" ")
    FILE_OUT.write("\n")          

FILE_OUT.close()

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

#To account for redundancies in NAME_ID
NAME_ID_FIXED=[]
for i in NAME_ID:
    if all(i[0]!=j[0] for j in NAME_ID_FIXED):
        NAME_ID_FIXED.append(i)
    
#FIX
FILE_IN=open("NETWORK/METABOLITES_ENDOGENOUS")
FILE=FILE_IN.read().splitlines()
FILE_IN.close()

FILE_OUT=open("NETWORK/METABOLITES_ENDOGENOUS_FIXED", "w")
for i in FILE:
    if i.split("$")[0].split("|")[1]!="None" and len(i.split("$")[0].split("|")[1])>3:
        for j in NAME_ID_FIXED:
            if i.split("$")[0].split("|")[0].lower()==j[0].lower():
                FILE_OUT.write(i.split("$")[0].split("|")[0]+"|"+j[1]+" "+"$"+" "+i.split("$")[1]+"\n")
    else:
        FILE_OUT.write(i+"\n")
FILE_OUT.close()   

