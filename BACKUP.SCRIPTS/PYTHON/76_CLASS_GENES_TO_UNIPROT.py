#CLASS GENES TO UNIPROT
class GENES_TO_UNIPROT_CLASS: #Takes list of genes and returns uniprots. Details of output in each definition
    "Implementation for the biomaRt (R) gene id to uniprot using the hsapiens_gene_ensembl"  
    
    #Load rpy2 needed libraries
    from rpy2 import robjects
    from rpy2.robjects import r
    
    #Declare r objects
    c=robjects.r["c"]
    #Load biomaRt
    r.library("biomaRt")
    #Choose biomaRt hsapiens
    ensembl=robjects.r["useMart"]("ensembl", dataset="hsapiens_gene_ensembl")
    
    def __init__(self, GENES):
        self.robjects=GENES_TO_UNIPROT_CLASS.robjects
        self.c=GENES_TO_UNIPROT_CLASS.c
        self.GENES=GENES
        self.ensembl=GENES_TO_UNIPROT_CLASS.ensembl
        #Vectorize gene list
        self.gene_vector=self.robjects.StrVector(self.GENES)
        #Call biomaRt
        RESULTS=self.robjects.r["getBM"](attributes=self.c("hgnc_symbol", "uniprot_swissprot_accession"), filters="hgnc_symbol",
                                         values=self.gene_vector, mart=self.ensembl)

        #Make RESULTS pythonic
        hgnc=list(RESULTS.rx2(1))
        uniprot=list(RESULTS.rx2(2))
        self.GENE_UNIPROT_PAIRS=[list(x) for x in zip(hgnc,uniprot)]
        #To filter out record that do not have a uniprot in file
        self.PRE_GENE_UNIPROT_PAIRS_FILTERED=filter(lambda x: len(x[1])>0, self.GENE_UNIPROT_PAIRS)
        
        #Further search against personal database
        LEFTOVER_GENE_SET=filter(lambda x: x not in zip(*self.PRE_GENE_UNIPROT_PAIRS_FILTERED)[0], self.GENES)
        import urllib2
        import pickle
        REVIEWED_FILE=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_UNIPROT_TO_GENES_REVIEWED.pi?attredirects=0&d=1")
        REVIEWED_DICT=pickle.load(REVIEWED_FILE)
        
        self.PLUS_GENE_UNIPROT_PAIRS=[]
        for gene in LEFTOVER_GENE_SET:
            for key,value in REVIEWED_DICT.iteritems():
                if gene in value:
                    self.PLUS_GENE_UNIPROT_PAIRS.append([gene, key])
                elif gene.lower() in value: #to account for lower case records, if already lower case and present earlier, it
                                            #would not get to this instance
                    self.PLUS_GENE_UNIPROT_PAIRS.append([gene, key])
                    
        #Add values found in personal database to biomart found values
        self.GENE_UNIPROT_PAIRS_FILTERED=self.PRE_GENE_UNIPROT_PAIRS_FILTERED+self.PLUS_GENE_UNIPROT_PAIRS
        
    def has_uniprot_genes (self): #Returns genes in query that have a uniprot match
        GENE_PAIRS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        self.found=[]
        for gene in self.GENES:
            if gene in GENE_PAIRS:
                self.found.append(gene)
        return self.found
    
    def get_uniprots (self): #Return list of all uniprots that correspond to uniprots in has_uniprot_genes
        self.UNIPROTS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[1]
        return self.UNIPROTS
    
    def no_uniprots(self): #Return genes for which no uniprot was found #CHECK BELOW FOR REASONS
        GENE_PAIRS=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        self.not_found=[]
        for gene in self.GENES:
            if gene not in GENE_PAIRS:
                self.not_found.append(gene)
        return self.not_found
    
    def not_on_database(self): #Return genes for which a record is not present in the ensembl database (hgnc) or personal
                                #database, double check them
        self.not_present=[]
        for gene in self.GENES:
            if gene not in zip(*self.GENE_UNIPROT_PAIRS)[0] and gene not in zip(*self.PLUS_GENE_UNIPROT_PAIRS)[0]:
                self.not_present.append(gene)
        return self.not_present
        
    def empty_gene(self): #Return genes that are in database but have no uniprot record - BE WARY OF RESULTS, relies on
                            #output format biomaRt, depending on size of input it may give "NA" or leave result empty giving the
                            #false impression that it is in the biomaRt database when in fact it does not exist
        self.empty=[]
        GENE_PAIRS_UNFILTERED=zip(*self.GENE_UNIPROT_PAIRS)[0]
        GENE_PAIRS_FILTERED=zip(*self.GENE_UNIPROT_PAIRS_FILTERED)[0]
        for gene in self.GENES:
            if gene in GENE_PAIRS_UNFILTERED and gene not in GENE_PAIRS_FILTERED:
                self.empty.append(gene)
        return self.empty

SEPARATOR=[]
#TEST RUN
import subprocess
LIST=subprocess.check_output("ls", 
                             cwd="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/geneStats").splitlines()

GENES=[]
for record in LIST:
    GENES.append(record.split(".")[0])
print GENES[:3]

#apply to class
print "number of genes in input:", len(GENES)
WANTED_GENES=GENES_TO_UNIPROT_CLASS(GENES)
print "genes that were found", len(WANTED_GENES.has_uniprot_genes()), WANTED_GENES.has_uniprot_genes()
print "genes that were not found", len(WANTED_GENES.no_uniprots()), WANTED_GENES.no_uniprots()
print "uniprots", len(WANTED_GENES.get_uniprots()), WANTED_GENES.get_uniprots()
print "not on database", len(WANTED_GENES.not_on_database()),WANTED_GENES.not_on_database()
print "on database but no known uniprot:", len(WANTED_GENES.empty_gene()), WANTED_GENES.empty_gene()
print "NA" in WANTED_GENES.get_uniprots()