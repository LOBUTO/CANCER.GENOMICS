#UNIPROT_CHAIN_LIMITS
#Returns list of mature start and end residue numbers given a uniprot
#If limits are unknown as in the following cases:
#    No Chain feature
#    ? chain number
#    < or > chain number
#The total length of the uniprot will be returned

def UNIPROT_CHAIN_LIMITS(UNIPROT_ID): #Given a uniprot, it returns the limits of the mature protein numbering
    import urllib2
    from Bio import SwissProt
    
    PAGE=urllib2.urlopen("http://www.uniprot.org/uniprot/%s.txt"%UNIPROT_ID)

    PARSED_PAGE=SwissProt.parse(PAGE)
    for record in PARSED_PAGE:
        CHAIN_VALUES=[]
        for feature in record.features:
            if feature[0]=="CHAIN":
                CHAIN_VALUES=CHAIN_VALUES+[str(feature[1]), str(feature[2])]

        if any(X.isdigit()==False for X in CHAIN_VALUES) or not CHAIN_VALUES:
            CHAIN_START=1
            CHAIN_END=record.sequence_length
        else:
            CHAIN_START=min(int(X) for X in CHAIN_VALUES)
            CHAIN_END=max(int(X) for X in CHAIN_VALUES)                
    
    return[CHAIN_START, CHAIN_END]  

A=UNIPROT_CHAIN_LIMITS("P00784")
print A, type(A)
