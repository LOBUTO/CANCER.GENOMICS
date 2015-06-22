#LEARN KEGG API

from urllib2 import urlopen
import os

#For pathway (hsa - human) have to look for both cpd(compound) and gl(glycan) 

#Get all human pathways and format to get "map" code
KEGG=urlopen("http://rest.kegg.jp/list/pathway/hsa").read().splitlines()
KEGG_MAP=["map"+X.split("\t")[0].split(":")[1].strip("hsa") for X in KEGG]

#Using "map" codes call cpd and gl 
GL=set([])
CPD=set([])
for record in KEGG_MAP:
    KEGG_GL_RECORD=urlopen("http://rest.kegg.jp/link/gl/%s"%record).read().splitlines()
    if len(KEGG_GL_RECORD)>1: #To account for empty records
        GLS=set([X.split("\t")[1].split(":")[1].strip() for X in KEGG_GL_RECORD])
        GL=GL | GLS
    
    KEGG_CPD_RECORD=urlopen("http://rest.kegg.jp/link/cpd/%s"%record).read().splitlines()
    if len(KEGG_CPD_RECORD)>1: #To account for empty records
        CPDS=set([Y.split("\t")[1].split(":")[1].strip() for Y in KEGG_CPD_RECORD])
        CPD=CPD | CPDS

print len(GL), GL
print len(CPD), CPD
