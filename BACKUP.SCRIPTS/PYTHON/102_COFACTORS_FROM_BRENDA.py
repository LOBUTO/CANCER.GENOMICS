#GET COFACTORS FOR EACH ENZYME (UNIPROT) USING BRENDA (SOAPpy) 
#USED SEVERAL TIMES TO AVOID DATABASE ERROR (Check notes in 8/21/13)
    # Need to add an all exception for this error later on

from FUNCTIONS import unique
from SOAPpy import SOAPProxy, Errors
import pickle, time, socket

#First get all EC numbers listed in BRENDA
BRENDA_URL="http://www.brenda-enzymes.org/soap2/brenda_server.php"
CLIENT=SOAPProxy(BRENDA_URL)

EC_NUMBERS=CLIENT.getEcNumbersFromEcNumber().split("!")
print EC_NUMBERS

#Second, get all uniprots per EC
UNIPROT_COFACTOR_IONS_DICT={}
for EC in EC_NUMBERS[5609:]: #if stuck at 4.1.99.18, start at [5236:]
    print EC    
    #Get uniprot list per EC
    while True:
        try:
            BRENDA_URL="http://www.brenda-enzymes.org/soap2/brenda_server.php"
            CLIENT=SOAPProxy(BRENDA_URL)
            UNIPROT_STRING=CLIENT.getSequence("ecNumber*%s#organism*Homo sapiens"%EC).split("!")
            if UNIPROT_STRING==[""]: UNIPROT_STRING=[]
            UNIPROTS=[]
            if (not UNIPROT_STRING)==False:
                UNIPROTS=[X.split("#")[3].split("*")[1] for X in UNIPROT_STRING]
                UNIPROTS=filter(lambda x: x!="more", UNIPROTS)
                print UNIPROTS
            break
        except Errors.HTTPError, e:
            if "502" in str(e):
                time.sleep(2)
                continue
            else:
                raise
        except socket.error, t:
            if "60" in str(t):
                time.sleep(60)
                continue
            else:
                raise
            
    #Get Cofactors per EC
    while True:
        try:
            BRENDA_URL="http://www.brenda-enzymes.org/soap2/brenda_server.php"
            CLIENT=SOAPProxy(BRENDA_URL)
            COFACTOR_STRING=CLIENT.getCofactor("ecNumber*%s"%EC).split("!") #ADD HOMO_SAPIENS FOR HUMAN ONLY!!!!!!!!
            if COFACTOR_STRING==[""]: COFACTOR_STRING=[]
            COFACTORS=[]
            if (not COFACTOR_STRING)==False:
                COFACTORS=list(set([X.split("#")[1].split("*")[1] for X in COFACTOR_STRING]))
                COFACTORS=filter(lambda y: y!="more", COFACTORS)
                print COFACTORS
            break
        except Errors.HTTPError, e:
            if "502" in str(e):
                time.sleep(2)
                continue
            else:
                raise
        except socket.error, t:
            if "60" in str(t):
                time.sleep(60)
                continue
            else:
                raise
            
    #Get Metals and ions
    while True:
        try:
            BRENDA_URL="http://www.brenda-enzymes.org/soap2/brenda_server.php"
            CLIENT=SOAPProxy(BRENDA_URL)
            METAL_ION_STRING=CLIENT.getMetalsIons("ecNumber*%s"%EC).split("!") #ADD HOMO_SAPIENS FOR HUMAN ONLY!!!!!!!!
            if METAL_ION_STRING==[""]: METAL_ION_STRING=[]
            METALS_IONS=[]
            if (not METAL_ION_STRING)==False:
                METALS_IONS=list(set([X.split("#")[1].split("*")[1] for X in METAL_ION_STRING]))
                METALS_IONS=filter(lambda z: z!="more", METALS_IONS)
                print METALS_IONS
            break
        except Errors.HTTPError, e:
            if "502" in str(e):
                time.sleep(2)
                continue
            else:
                raise
        except socket.error, t:
            if "60" in str(t):
                time.sleep(60)
                continue
            else:
                raise
                        
    #Mix Side products (Cofactors, metals and ions)
    SIDE_PRODUCTS=COFACTORS+METALS_IONS
    
    #Create dict per UNIPROT
    if len(UNIPROTS)>0 and len(SIDE_PRODUCTS)>0:
        RECORD_DICT=dict((X,SIDE_PRODUCTS) for X in UNIPROTS)
        UNIPROT_COFACTOR_IONS_DICT.update(RECORD_DICT)
    
    print UNIPROT_COFACTOR_IONS_DICT
    
    #Pickle at every round in case process runs out due to unsolvable error
    PICKLE_OUT=open("DATABASES/OBJECTS/DICT_UNIPROT_TO_SIDE_PRODUCT_ALL_PART9.pi", "w")
    pickle.dump(UNIPROT_COFACTOR_IONS_DICT, PICKLE_OUT)
    PICKLE_OUT.close()