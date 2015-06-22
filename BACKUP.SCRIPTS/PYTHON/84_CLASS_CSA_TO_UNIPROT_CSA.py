#CLASS TO PARSE CSA SQL TO UNIPROT-CS DICTIONARY OBJECT 7/27/13

class CSA_TO_UNIPROT_CSA_CLASS: #Dictionary Object - Make sure server is running an you have imported SQL dump file into
                                        #a database in MYSQL 
    #Load functions
    import MySQLdb as MQ
    
    def __init__ (self, SQL_DATABASE, HOST, USER, PASSWORD, LOG_FOLDER):
        from FUNCTIONS import unique, PDB_RES_TO_UNIPROT_RES        
        self.unique=unique
        self.PDB_RES_TO_UNIPROT_RES=PDB_RES_TO_UNIPROT_RES
        self.SQL_DATABASE=SQL_DATABASE
        self.HOST=HOST
        self.USER=USER
        self.PASSWORD=PASSWORD
        
        #Create connection object that represents database
        conn_obj=CSA_TO_UNIPROT_CSA_CLASS.MQ.connect(host=self.HOST,user=self.USER, passwd=self.PASSWORD, db=self.SQL_DATABASE)
        cursor=conn_obj.cursor()
        
        #Query Catalytic Residues table for PDB_ID, chain, PDB_residue_number, uniprot_residue_number, residue type 
        cursor.execute("SELECT PDBID, CHAIN, NUMBER TYPE FROM CSA_CATALYTICRESIDUES")
        PRE_RESULTS=[list(x) for x in cursor.fetchall()]
        PRE_RESULTS=self.unique(PRE_RESULTS)
        print PRE_RESULTS
        
        #Format for PDB_RES_TO_UNIPROT_RES function
        PDB_DICT=dict((x[0].upper()+"|"+x[1],[]) for x in PRE_RESULTS)
        print PDB_DICT
        for record in PRE_RESULTS:
            PDB_DICT[record[0].upper()+"|"+record[1].upper()]=PDB_DICT[record[0].upper()+"|"+record[1].upper()]+[int(record[2])]
        print PDB_DICT
        
        #Call function
        self.PDB_UNIPROT_RESULTS=self.PDB_RES_TO_UNIPROT_RES(PDB_DICT, LOG_FOLDER)

        #Find those that have UNIPROT-RES answer (match)
        self.WITH_ALL=dict((X, self.PDB_UNIPROT_RESULTS[X]) for X in self.PDB_UNIPROT_RESULTS.iterkeys() 
                           if self.PDB_UNIPROT_RESULTS[X]!="None" and len(self.PDB_UNIPROT_RESULTS[X])!=0) 
            
    def WITH_MATCH_ALL (self): #Returns dictionary in the form {"PDB|CHAIN=UNIPROT":[csa_residues]} if PDB had a match
        return self.WITH_ALL
        
    def WITH_MATCH_UNIPROT (self): #Returns uniprot dictionary of those PDB that had a match in the form {"UNIPROT""[csa_residues]}
        self.WITH_ALL_UNIPROT=dict((X.split("=")[1], self.WITH_ALL[X]) for X in self.WITH_ALL.iterkeys())
        
        return self.WITH_ALL_UNIPROT
        
    def WITHOUT_MATCH (self): #Returns dictionary in the form {"PDB|X":"None"} for those PDB that did not match a uniprot
        self.WITH_NONE=dict((X, self.PDB_UNIPROT_RESULTS[X]) for X in self.PDB_UNIPROT_RESULTS.iterkeys() 
                           if self.PDB_UNIPROT_RESULTS[X]=="None" or len(self.PDB_UNIPROT_RESULTS[X])==0)          
        return self.WITH_NONE
    
#TEST
TEST=CSA_TO_UNIPROT_CSA_CLASS("TESTING", "localhost", "root", "mysql", "LOGS")
print "all",TEST.WITH_MATCH_ALL()
print "uniprot", TEST.WITH_MATCH_UNIPROT()
print "without", TEST.WITHOUT_MATCH()
print  TEST.WITH_MATCH_UNIPROT()["    "]
print TEST.WITH_MATCH_UNIPROT().has_key("")
print TEST.WITH_MATCH_UNIPROT().has_key("P22643")