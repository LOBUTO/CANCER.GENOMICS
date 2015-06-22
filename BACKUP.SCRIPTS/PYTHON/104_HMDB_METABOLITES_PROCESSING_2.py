#HMDB METABOLITES PROCESSING 2 - HELPER
#08/15/13
#Process DRUGBANK.xml to get tags

import xml.etree.ElementTree as ET

#Process file
ROOT=ET.parse("DATABASES/drugbank.xml").getroot()

#Get names
NAMES=[X[1].text for X in ROOT[:]]

#Write to output
FILE_OUT=open("DATABASES/DRUGBANK_NAMES", "w")
for record in NAMES:
    FILE_OUT.write(record+"\n")
FILE_OUT.close()

#Processed output afterwards to delete some extra lines!!!!!!!!!!!!