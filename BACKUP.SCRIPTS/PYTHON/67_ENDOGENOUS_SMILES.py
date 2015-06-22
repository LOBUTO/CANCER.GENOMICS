#GET COMBINED smiles FILE FOR ALL ENDOGENOUS METABOLITES (FROM METABOLITES_ENDOGENOUS_FIXED, FILE 53.py)

import xml.etree.ElementTree as ET
import subprocess

#Get metabolite_names|PDB_ID
FILE_IN=open("NETWORK/METABOLITES_ENDOGENOUS_FIXED").read().splitlines()
END_MET=[]
for line in FILE_IN:
    END_MET.append(line.split("$")[0])
print len(END_MET)

#Get original smiles files dictionary
XML_FILES=subprocess.check_output("ls", cwd="DATABASES/hmdb_metabolites").splitlines()

SMILES_DICT={}
for xml in XML_FILES:
    ROOT=ET.parse("DATABASES/hmdb_metabolites/"+xml).getroot()
    for origin in ROOT.iter("origin"):
        if origin.text=="Endogenous" and str(ROOT.find("smiles").text)!="None":
            SMILES_DICT[ROOT.find("name").text]=ROOT.find("smiles").text

#Get terms
FILE_OUT=open("/Users/jzamalloa/Desktop/FOLDER/LOGS/END_SMILES", "w")
END_SMILES=[]
for met in END_MET:
    if SMILES_DICT.has_key(met.split("|")[0]):
        END_SMILES.append([SMILES_DICT[met.split("|")[0]], met.split("|")[0]])
        FILE_OUT.write(SMILES_DICT[met.split("|")[0]]+"    "+met.split("|")[0]+"\n")
    else:
        print met
FILE_OUT.close()
print len(END_SMILES)