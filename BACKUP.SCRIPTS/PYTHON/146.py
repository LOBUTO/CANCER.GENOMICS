#032414
#Get enyzme/product table from kegg

from urllib2 import urlopen
import itertools

#Get all human (hsa) KEGG genes IDs and their respective Hugo symbols
URL_IN1=urlopen("http://rest.kegg.jp/list/hsa").read().splitlines()
DICT_GENE=dict((x.split("\t")[0].split(":")[1],
                x.split("\t")[1].split(",")[0].split(";")[0].strip()) for x in URL_IN1)

print len(DICT_GENE) #30687 human enzymes, including those with uncharacterized names

#Get all compounds CPD:[Names]
URL_IN2=urlopen("http://rest.kegg.jp/list/cpd/").read().splitlines()
DICT_CPD=dict((y.split("\t")[0].split(":")[1],
               [z.strip() for z in y.split("\t")[1].split(";")]) for y in URL_IN2)

#Get all glycans and add them to the compounds dictionary
URL_IN3=urlopen("http://rest.kegg.jp/list/gl/").read().splitlines()
DICT_GL=dict((y.split("\t")[0].split(":")[1],
               [z.strip() for z in y.split("\t")[1].split(";")]) for y in URL_IN3)
DICT_CPD.update(DICT_GL)

#Get all enzymes - Filtering for those that are obsolete
URL_IN4=urlopen("http://rest.kegg.jp/list/enzyme").read().splitlines()
ENZYMES=[y.split("\t")[0] for y in URL_IN4 if y.split("\t")[1][:14] !="Transferred to"] #5584 Enzyme codes

#Write to table
FILE_OUT1=open("DATABASES/KEGG/032414_ENZYME_PRODUCT_TEST","w") #PRODUCT_TEST tested run on 061214, changed name to 061214
FILE_OUT1.write("Enzyme"+"\t"+"KEGG_ID"+"\t"+"Product")

#Process and build table
for enzyme in ENZYMES:
    print enzyme
    TEMP_URL_IN1=urlopen("http://rest.kegg.jp/link/genes/"+enzyme).read().splitlines()
    
    #Check that enzyme page has an actual annotated gene before continuing
    if TEMP_URL_IN1[0]!="":
        
        #Get hsa KEGG genes ids
        HSAS=[z.split("\t")[1].split(":")[1] for z in TEMP_URL_IN1 if z.split("\t")[1].split(":")[0]=="hsa"]
        
        #Continue only if it has human genes (hsa)
        if len(HSAS)>0:
            
            #Translate the gene KEGG id to their Hugo identifier
            HSAS_TRANSLATE=[DICT_GENE[Z] for Z in HSAS]
            
            #Get reactions for this enzyme
            TEMP_URL_IN2=urlopen("http://rest.kegg.jp/link/rn/"+enzyme, timeout=3000).read().splitlines()
            
            #Check that enzyme page has an actual annotated reaction before continuing
            if TEMP_URL_IN2[0]!="":
            
                REACTION=[w.split("\t")[1] for w in TEMP_URL_IN2]
                
                for reaction in REACTION:
                    print reaction
                    
                    #Scrape the url for each reaction to get the products compound codes
                    TEMP_URL_IN3=urlopen("http://rest.kegg.jp/get/"+reaction)
                    PRODUCTS=[x.strip("EQUATION").strip().split("<=>")[1].strip().split(" + ") for x in TEMP_URL_IN3 if x[0:8]=="EQUATION"][0]
                    
                    ##Clean stochiometric numbers from compounds - LOOK IN THE FUTURE TO SEE HOW WE CAN USE THIS*
                    PRODUCTS=list(itertools.chain(*[x.split() for x in PRODUCTS]))
                    PRODUCTS=filter(lambda y: len(y)>5, PRODUCTS)
                    PRODUCTS=[x.replace("(n-1)","").replace("(n)","").replace("(n+1)","").replace("(n-2)","").replace("(n+2)","").replace("(m)","").replace("(n+m)","") 
                              for x in PRODUCTS]
                    
                    #Write to table for each synonym of compound and each gene in enzyme category
                    for cpd in PRODUCTS:
                        print cpd
                        
                        cpd_names=DICT_CPD[cpd]
                        
                        for name in cpd_names:
                            for gene in HSAS_TRANSLATE:
                                FILE_OUT1.write("\n"+gene+"\t"+cpd+"\t"+name)

FILE_OUT1.close()
                
#NEEDS POST-PROCESSING!!
#Such as unique()