#UNIPROT XML PROCESSING SAMPLE
#01/18/14


import xml.etree.cElementTree as ET

ROOT=ET.parse("DATABASES/UNIPROT/uniprot-%28organism%3A%22Homo+sapiens+%5B9606%5D%22%29+AND+reviewed%3Ayes.xml").getroot()

FILE_OUT1=open("DATABASES/UNIPROT/GENE_LOCATIONS", "w")
FILE_OUT1.write("GENE"+"\t"+"LOCATION")

for record in ROOT:
    
    #Get GENE name - Some of the records don't have a gene assignemt
    GENE_NAME=record.find("{http://uniprot.org/uniprot}gene/{http://uniprot.org/uniprot}name") 
    if GENE_NAME!=None:
        GENE_NAME=GENE_NAME.text.strip()
    
        LOCATIONS=record.findall("{http://uniprot.org/uniprot}comment/{http://uniprot.org/uniprot}subcellularLocation/{http://uniprot.org/uniprot}location")
        LOCATIONS=[X.text for X in LOCATIONS]
    
        for location in LOCATIONS:
            FILE_OUT1.write("\n"+GENE_NAME+"\t"+location)
        
FILE_OUT1.close()
        
    