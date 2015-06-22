#ADDITIONAL PDBS FOR PDB_UNIPROT OBJECT
from prody import *
from FUNCTIONS import GET_PDB
import subprocess
import os
import pickle

#List of PDBs that we need a uniprot for:
PDBS= ['1a65', '1b65','1bqc', '1c4x', '1d2t', '1dio', '1ex1', '1ga8', '1h4g', '1h54', '1i19', '1im5', '1it4', '1jfl', '1jkm', '1jnr', '1jrp', '1kzh', '1l1r', '1m21', '1mpx', '1mqw', '1n2t', '1nf9', '1nml', '1okg', '1p4n', '1pj5', '1pjb', '1pyl', '1qaz', '1qgn', '1r4f', '1r76', '1rhc', '1ro7', '1tz3', '1ujn', '1uk7', '1un1', '1v0y', '1x7d', '1xny', '1y9m', '1zoi', '2ag0', '2bsx', '2c7v', '2dbt', '2dw7', '2nlr', '2pda', '2qf7']

#for pdb in PDBS:
#    GET_PDB(pdb, "PDB_FILES")
"""
FILES=subprocess.check_output("ls", cwd="PDB_FILES").splitlines()
PDB_FILES=[x for x in FILES if "ent" in x]
print PDB_FILES
print len(PDB_FILES)

os.chdir("PDB_FILES")
for item in PDB_FILES:
    subprocess.Popen(["-c", "mv %s %s"%(item, item[3:7])],stdout=subprocess.PIPE, shell=True)
"""
SEPARATOR=[]
#Use urllib2,urllib on http://www.uniprot.org/mapping/ to map pdb->uniprot as tab table
#Save list as an object in the format: [[pdb, [uniprot1, uniprot2,..],...]

def PDB_TO_UNIPROT(PDB_LIST): #DICT of format {pdb:[unipro1, uniprot2,...],...}
                                #Uses http://www.uniprot.org/mapping/ server to retrieve query results as dictionary.
                                #WARNING, Only use when you are sure of the validity of your pdb list, a better choice that is
                                #curated is the class PDB_TO_UNIPROT_CLASS, use this as supplementary to the class. This function
                                #may not tell you if pdbs have been replaced by other while querying by server, use with caution 
    import urllib, urllib2

    url="http://www.uniprot.org/mapping/"

    #Set up query - It has to be as a single string of elements
    QUERY=" ".join(PDB_LIST)

    #Set up values, attributes and filters 
    params={'from':"PDB_ID", 'to':"ACC", 'format':"tab", 'query':QUERY}

    #Send job
    DATA=urllib.urlencode(params) #Encode parameters and values for form
    REQUEST=urllib2.Request(url, DATA) #Form
    RESPONSE=urllib2.urlopen(REQUEST) #Get results
    PAGE=RESPONSE.read(200000) #Process Results object into readable form

    #Process Page response
    PDB_UNIPROT_LIST=[X.split("\t") for X in PAGE.splitlines()[1:]]
    
    #Build dictionary
    PDB_UNIPROT_DICT=dict((x.upper(), []) for x in zip(*PDB_UNIPROT_LIST)[0])
    for record in PDB_UNIPROT_LIST:
        PDB_UNIPROT_DICT[record[0].upper()]=PDB_UNIPROT_DICT[record[0].upper()]+[record[1]]
    
    return PDB_UNIPROT_DICT

TEST=PDB_TO_UNIPROT(PDBS)    
print TEST

PICKLE_OUT=open("DATABASES/OBJECTS/DICT_PDB_TO_UNIPROT_ALL_SUP.pi", "w")
pickle.dump(TEST,PICKLE_OUT)
