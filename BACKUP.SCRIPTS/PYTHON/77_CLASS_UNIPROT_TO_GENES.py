#UNIPROT_TO_GENES_CLASS
#Returns uniprot results given a list of genes
#08/06/13
#Given a list of UNIPROTs returns main results in uniprot_gene module as dict {UNIPROT:[GENE1, GENE2,...]}
#Non matches are shown as a list in the no-genes module

class UNIPROT_TO_GENES_CLASS: #Takes in list UNIPROT IDs and returns hgnc gene symbols 
    "loading class of uniprot-to-gene dictionary"
    
    import urllib2
    import pickle
    FILE_IN=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_UNIPROT_TO_GENES.pi?attredirects=0&d=1")
    DICT=pickle.load(FILE_IN) #Loads the uniprot_to_genes dictionary
    
    def __init__(self, UNIPROT):
        self.uniprot=UNIPROT
        self.dict=UNIPROT_TO_GENES_CLASS.DICT
    
    def uniprot_gene (self): #Returns dictionary of the form {UNIPROT:[GENE1, GENE2....]} for uniprots that have a matching gene
        self.UNIPROT_DICT={}
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                self.UNIPROT_DICT[record]=self.dict[record]
        return self.UNIPROT_DICT
    
    def has_genes_uniprots(self): #Returns lists of uniprots that have genes
        TRUE_UNIPROTS=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                TRUE_UNIPROTS.append(record)
        return TRUE_UNIPROTS
    
    def get_genes(self): #Returns list of all the genes that are present in given uniprots that have keys and are in database
        ALL_AVAILABLE_GENES=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])!=0:
                ALL_AVAILABLE_GENES.extend(self.dict[record])
        return ALL_AVAILABLE_GENES
        
    def no_genes(self): #Returns lists of uniprots for which genes cannot be retrieved, CHECK WITH BELOW FUNCTIONS FOR WHY
        NO_UNIPROT_GENES=[]
        for record in self.uniprot:
            if self.dict.get(record)==None or len(self.dict[record])==0:
                NO_UNIPROT_GENES.append(record)
        return NO_UNIPROT_GENES
    
    def not_on_file(self): #Returns list of uniprots that are not present in list - IMPORTANT TO ADD THEM TO FILE
        ABSENT_UNIPROTS=[]
        for record in self.uniprot:
            if self.dict.get(record)==None:
                ABSENT_UNIPROTS.append(record)
        return ABSENT_UNIPROTS
    
    def empty(self): #Returns list of uniprots that are on file but have no annotated genes - #DEPRECATED
        self.EMPTY=[]
        for record in self.uniprot:
            if self.dict.has_key(record) and len(self.dict[record])==0:
                self.EMPTY.append(record)
        return self.EMPTY

SEPARATOR=[]

#Load list of endogenous uniprots
FILE_IN1=open("NETWORK/ENDOGENOUS_UNIPROT").read().splitlines()
print FILE_IN1[:4]
print len(FILE_IN1)

UNIPROTS_CLASS=UNIPROT_TO_GENES_CLASS(FILE_IN1)
print "have genes:", len(UNIPROTS_CLASS.has_genes_uniprots()) 
print "don't have genes", len(UNIPROTS_CLASS.no_genes()), UNIPROTS_CLASS.no_genes()
for uniprot in UNIPROTS_CLASS.no_genes(): print uniprot
print "not on list", UNIPROTS_CLASS.not_on_file()
print "on file but not available genes", UNIPROTS_CLASS.empty()
print "RESULTS", UNIPROTS_CLASS.uniprot_gene()