#WORKING WITH PRINCETON GENE ONTOLOGY ENRICHMENT SITE @ http://go.princeton.edu/cgi-bin/GOTermFinder
#Function to determine GO enrichment from gene list
#BUILT 06/4/13 - IMPORTANT - If site format changes then program may have unexpected and unrealible consequences, double check script for
#accuracy from time to time

#Go to site
def GO_ENRICHMENT(GENE_LIST_FILE, BACKGROUND_FILE, FILE_OUT,GENE_ONTOLOGY="Process" ,p_cutoff=0.01): #GENE_LIST_FILE can be a single one 
                                                                #gene-per-line file or for batch jobs can be zipped files of this format, 
                                                                #BACKGROUN_FILE is a single gene-per-line file, FILE_OUT is where results
                                                                #will be stored, GENE_ONTOLOGY can either be: "Process", "Function" or 
                                                                #"Component", p_cutoff is the p-value cutoff for significant shared GO terms
                                                                #search                                                                
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
                                                                                    
    #Go to site
    DRIVER=webdriver.Firefox()
    DRIVER.set_window_size(600,400)    
    DRIVER.get("http://go.princeton.edu/cgi-bin/GOTermFinder")
    assert "Gene Ontology" in DRIVER.title
    
    #Introduce list of genes
    UPLOAD_GENE=DRIVER.find_element_by_name("uploadFile") #could have chosen to find it by other means, chose by_name, as long as finding
                                                            #parameter was within the "input bracket" in html obtained by page_source
    UPLOAD_GENE.send_keys(GENE_LIST_FILE)

    #Choose and ontology:
    GO=DRIVER.find_elements_by_tag_name("label") #In this instance we choose "label" because it is a checkbox, similar to select
    for j in GO:
        if j.text.lower()==GENE_ONTOLOGY.lower():j.click()
    
    #Choose annotation
    ANNOTATIONS=DRIVER.find_elements_by_tag_name("option") #again could have chosen other, but since "option" is a parameter within each 
                                                            #element of drop down menu chose it instead 
    ORGANISMS=[]
    X=1
    for organism in ANNOTATIONS: 
        ORGANISMS.append(str(X)+"#"+organism.text)
        X=X+1
    for i in ORGANISMS: print i.split("#")[0], i.split("#")[1]
    
    USER_SELECTION2=raw_input("Choose annotation by entering number key and enter:")
    
    SELECTION_ORGANISM=filter(lambda x:x.split("#")[0]==USER_SELECTION2, ORGANISMS)
    SELECTION_ORGANISM=SELECTION_ORGANISM[0].split("#")[1]
    print SELECTION_ORGANISM

    for option in ANNOTATIONS:
        if option.text==SELECTION_ORGANISM: option.click()
    
    #Upload background
    BACKGROUND=DRIVER.find_element_by_name("uploadBackground")
    BACKGROUND.send_keys(BACKGROUND_FILE)
    
    #Set p-value cut-off - 0.01 is used as default
    P_VALUE=DRIVER.find_element_by_name("pValueCutoff")
    P_VALUE.send_keys(str(p_cutoff))
    
    #SUBMIT
    SUBMIT=DRIVER.find_element_by_name("searchGOTermsButton")
    SUBMIT.click()
    
    #Check for completion and mistakes in input
    Z=0
    while Z==0:
        URL=DRIVER.current_url
        CURRENT_PAGE=[]
        for s in DRIVER.page_source.splitlines():
            CURRENT_PAGE.append(s)
        
        if "JobWatch" in URL: 
            if len(CURRENT_PAGE)==16:#WARNING - This is extremely dependent on the output format
                                     #so test with visuals every once in a while
                Z=2        
            elif len(CURRENT_PAGE)==17:
                Z=0
                URL=DRIVER.current_url
        else:
            Z=1

    #If no problems with run
    print Z
    if Z==1:
        #GET OUTPUT:
        WEB_RESULT_LINK=DRIVER.find_element_by_link_text("Tab-delimited")#This time we find by link, if you look at the HTTP it is difficult to classify
                                                                        #the link by a parameter, that's how you know it could be a link
        WEB_RESULT_LINK.click()
        
        #SAVE
        WEB_RESULT=[]
        for string in DRIVER.page_source.splitlines():
            WEB_RESULT.append(string)
        WEB_RESULT=WEB_RESULT[1:-1]
        
        OUTPUT=open(FILE_OUT,"w")
        for line in WEB_RESULT:
                OUTPUT.write(line+"\n")
                
        print "File saved to %s"%FILE_OUT
    else:
        print "No known genes found. Check that genes in list are annotated with GO terms and that GO association (organism) is correct"
    
    #Close browser
    DRIVER.quit() 
        

DUMMY=[]
#TEST
GO_ENRICHMENT("/Users/jzamalloa/Desktop/GO-TermFinder-0.86/examples/genes.txt", 
              "/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/test",
              "GO/LOGS/TEST", "Function",0.05)


