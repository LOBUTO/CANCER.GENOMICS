#Learning urllib

import urllib,urllib2

url=urllib2.urlopen("http://www.uniprot.org/uniprot/F1PBL3.fasta")
URLA=[]
for i in url: URLA.append(i)
print URLA
URLB=""
for j in URLA[1:]:
    URLB=URLB+j[:-1]
print URLB
    


def UNIPROT_SEQ(UNIPROT): #STRING - Feed UNIPROT ID and returns FASTA sequence
    SEQ=""
    url=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.fasta"%UNIPROT).readlines()
    for i in url[1:]:
        SEQ=SEQ+i[:-1]
    return SEQ

UNI="F1PBL3"

print 7
CALL=UNIPROT_SEQ(UNI)
print CALL
