#LEARNING PRODY
from prody import *
"""
atom, header=parsePDB("/Users/jzamalloa/Downloads/1A4Y.pdb" ,header=True)

print list(header)
print list(header["polymers"])
for row in list(header["polymers"]):
    print row

print header["A"]
print header["B"]
print header["identifier"]
print header["classification"]
print 7
for i in list(header):
    print header[i]
"""    
print "PARSING HEADERS"  
#PARSING HEADERS
POLYMER=parsePDBHeader("PDB_FILES/ALL_PDB/1bcc", "polymers")

for polymers in POLYMER: #BUILD LISTS TO QUERY BY .chid for CSA access
    print polymers.chid, polymers.name
    print polymers.dbrefs[0].database=="UniProt", polymers.dbrefs[0].accession #Keep in mind that this work if database we are looking at
                                                                    #is uniprot, so you have to match it before getting it
                                                                    

