#CLASS FOR GETTING UNIPROT GIVEN A PDB LIST

class PDB_TO_UNIPROT_CLASS: #Given a pdb list it produces uniprot results
    import urllib2
    import pickle
    
    #Load PDB_TO_UNIPROT personal dictionarys
    DICT_FILE1=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_TO_UNIPROT_ALL.pi?attredirects=0&d=1")
    DICT=pickle.load(DICT_FILE1)
    
    #THIS IS USED FOR FUTURE UPDATES
    DICT_FILE2=urllib2.urlopen("https://sites.google.com/site/jdatabases/home/uniprot/DICT_PDB_TO_UNIPROT_ALL_SUP1.pi?attredirects=0&d=1")
    DICT.update(pickle.load(DICT_FILE2))
    ################################
    
    def __init__(self,PDB_LIST):
        self.DICT=PDB_TO_UNIPROT_CLASS.DICT
        self.PDB_LIST=PDB_LIST
        
    def has_uniprot_pdbs (self): #Displays pdbs that have a uniprot in record
        self.with_uniprot=filter(lambda x: self.DICT.has_key(x.upper())==True, self.PDB_LIST)
        return self.with_uniprot    
        
    def no_uniprot_pdbs (self): #Displays pdbs that do not have a uniprot in record
        self.without_uniprot=filter(lambda x: self.DICT.has_key(x.upper())==False, self.PDB_LIST)
        return self.without_uniprot
    
    def uniprots (self): #Displays list of matching uniprots of pdb that have uniprots #KEEP IN MIND COMPLEXES WILL HAVE FORMAT
                            #*|*
        self.uniprot_list=["|".join(self.DICT[X.upper()]) for X in self.with_uniprot]
        return self.uniprot_list
    
    def uniprot_pdb(self): #Returns list of pair [pdb,uniprot] of those pds that have a uniprot on record -FORMAT AS ABOVE
        self.results={}
        for pdb in self.with_uniprot:
            self.results[pdb.upper()]="|".join(self.DICT[pdb.upper()])
        return self.results
        
SEPARATOR=[]
PDB_L=["10GS", "NOPP", "11GS","1bcc", "1a65"]
TEST=PDB_TO_UNIPROT_CLASS(PDB_L)
print TEST.has_uniprot_pdbs()
print TEST.no_uniprot_pdbs()
print TEST.uniprots()
print TEST.uniprot_pdb()