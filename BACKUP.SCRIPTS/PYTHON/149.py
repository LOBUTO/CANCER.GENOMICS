#Funciton to get KEGG cpd per pathway
#050314
#Remember to make into function!!

from urllib2 import urlopen

#Load human pathways from KEGG
URL_IN1=urlopen("http://rest.kegg.jp/list/pathway/hsa").read().splitlines()
PATHWAYS_DICT=dict((X.split()[0].split(":")[1].strip("hsa"), X.split("\t")[1]) for X in URL_IN1)


#Create file to write out to
#FILE_OUT1=open("DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND", "w")
FILE_OUT1=open("DATABASES/KEGG/060415_PATHWAY_TO_COMPOUND", "w")
FILE_OUT1.write("PATHWAY"+"\t"+"COMPOUND"+"\t"+"DESCRIPTION")

#To keep track
W=0

#Get compounds and glycans per pathway
for path in PATHWAYS_DICT.keys():
    
    #Get compounds first:
    URL_IN2=urlopen("http://rest.kegg.jp/link/cpd/map"+path).read().splitlines()
    
    #Check that there are actually compounds for pathway 
    if URL_IN2[0]!="":
        
        #Get compounds
        CPDS=[y.split()[1].split(":")[1] for y in URL_IN2]
        
        #Write to table
        for compound in CPDS:
            FILE_OUT1.write("\n"+path+"\t"+compound+"\t"+PATHWAYS_DICT[path])
    
    #Then get glycans:
    URL_IN3=urlopen("http://rest.kegg.jp/link/gl/map"+path).read().splitlines()
    
    #Check that there are actually glycans for pahtway
    if URL_IN3[0]!="":
        
        #Get compounds
        GLS=[z.split()[1].split(":")[1] for z in URL_IN3]
        
        #Write to table
        for glycans in GLS:
            FILE_OUT1.write("\n"+path+"\t"+glycans+"\t"+PATHWAYS_DICT[path])
    
    W=W+1
    print float(W)/len(PATHWAYS_DICT)
FILE_OUT1.close() 