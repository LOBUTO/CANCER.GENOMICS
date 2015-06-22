#LEARNING XML PARSING WITH ElementTree

import xml.etree.ElementTree as ET

#To parse the file so you could start iterating and using it as a tree (nested lists)
FILE=ET.parse("DATABASES/hmdb_metabolites/HMDB00001.xml")

#Very important to get the root so you have somewhere to start and iterate over all of the nested children
ROOT=FILE.getroot()
print ROOT
print ROOT.tag #Prints the tag, you could say this is the dictionary key of the level 
print ROOT.attrib #Prints the attrib which is a dictionary type, so if you had an "= ...." at the end of the level lines instead of < > then
                    #it would print is as a dictionary term 

#Get each element in the lower level of what you are iterating over 
for child in ROOT: 
    print child.tag, child.attrib
    
#Since it is an iterable tree, you could treat it as list of multiple level and go through it
print ROOT[4][2].tag

#To print the acutal value you use .text
print ROOT[4][1].text

print 7
#To iterate over it (including children and children of children) and find actual information (keys, values) use .iter("term")
for gene_name in ROOT.iter("gene_name"):
    print gene_name.text

#IMPORTANT - Haven't yet found a way to get specific children without iterating through rooted xml

#UPDATE
print ROOT.find("name").text #Finds the first subelement (from level being called) only, if want a list of all not just first use findall()
print type(ROOT.find("name").text)  