#Learning BRENDA SOAPpy

import string, socket, time
from SOAPpy import SOAPProxy, WSDL

while True:
    try:
        endpointURL = "http://www.brenda-enzymes.org/soap2/brenda_server.php"
        client = SOAPProxy(endpointURL)

        print 7
        #UNIPROT AC
        resultString = client.getSequence("ecNumber*1.1.1.1#organism*Homo sapiens").split("!")
        print resultString
        for record in resultString:
            print record.split("#")
        
        ACCESSIONS=[X.split("#")[3].split("*")[1] for X in resultString]
        print ACCESSIONS
        print 7
        #COFACTOR
        print client.getCofactor("ecNumber*4.1.99.18#organism*Homo sapiens")
        COFACTOR_STRING=client.getCofactor("ecNumber*4.1.99.18#organism*Homo sapiens").split("!")
        print COFACTOR_STRING, len(COFACTOR_STRING)
        
        COFACTORS=[X.split("#")[1].split("*")[1] for X in COFACTOR_STRING]
        print COFACTORS
        break
        
    except socket.error, t:
        if "60" in str(t):
            time.sleep(2)
            continue
        else:
            raise   