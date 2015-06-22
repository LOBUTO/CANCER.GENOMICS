#PDB_CHAIN_TO_UNIPROT CLASS - 8/3/13
#Dependencies:
#    PDB_TO_UNIPROT_CLASS
#    BLAST_PDB_AGAINST_CLASS
#Returns dictionary of {PDB|X:UNIPROT}
#If PDB|X does not exist it will return error 
#If there is absolutely no matching UNIPROT it will return "None"
class PDB_CHAIN_TO_UNIPROT_CLASS: #Given a LIST!!! of PDB|X returns Uniprots corresponding to it
    
    def __init__ (self, PDB_LIST):
        #Load modules
        import subprocess
        from FUNCTIONS import PDB_TO_UNIPROT_CLASS, BLAST_PDB_AGAINST_CLASS
        
        #Load list
        self.PDB_LIST=[X.upper() for X in PDB_LIST]
        
        #Check for pdbs using class
        CLASS1=PDB_TO_UNIPROT_CLASS([X.split("|")[0] for X in self.PDB_LIST])
        CLASS_HAS=CLASS1.has_uniprot_pdbs()
        CLASS_HASNT=CLASS1.no_uniprot_pdbs()

        #For those pdbs that have uniprots, targeted blast function 
        self.PDB_LIST_HAS=[X for X in self.PDB_LIST if X.split("|")[0] in CLASS_HAS]
        self.PDB_CHAIN_UNIPROT={} #Dictionary to store those PDB|X that do have a uniprot
        self.PDB_CHAIN_NO_UNIPROT={} #Dictionary to store those PDB|X that do not have a uniprot, value="None"
        
        for record in self.PDB_LIST_HAS:
            CLASS2=BLAST_PDB_AGAINST_CLASS(record, 60, 80, "LOGS")
            UNIPROT_CANDIDATES=CLASS1.uniprot_pdb()[record.split("|")[0]].split("|") #Obtained from PDB_TO_UNIPROT Module
            CLASS2_RESULT=CLASS2.TO_SUBJECT(UNIPROT_CANDIDATES)

            if CLASS2_RESULT!="None":
                self.PDB_CHAIN_UNIPROT[record]=CLASS2_RESULT
            else:
                self.PDB_CHAIN_NO_UNIPROT[record]="None"
                
        #For those pdbs that do not have uniprot total blast function
        self.PDB_LIST_HASNT=[X for X in self.PDB_LIST if X.split("|")[0] in CLASS_HASNT]
        
        for record in self.PDB_LIST_HASNT:
            CLASS3=BLAST_PDB_AGAINST_CLASS(record, 60, 80, "LOGS")
            CLASS3_RESULT=CLASS3.TO_ALL()
            
            if CLASS3_RESULT!="None":
                self.PDB_CHAIN_UNIPROT[record]=CLASS3_RESULT
            else:
                self.PDB_CHAIN_NO_UNIPROT[record]="None"
    
    def with_uniprots (self): #Returns dictionary of PDB|X:uniprot that do have single qualifier uniprot
        return self.PDB_CHAIN_UNIPROT
        
    def without_uniprots (self): #Returns dictionary of PDB|X:"None" that do not have a uniprot
        return self.PDB_CHAIN_NO_UNIPROT

#TEST CLASS
PDBS=["10GS|A", "11GS|B","1bcc|A", "1a65|a", "9PAP|A"]
TEST=PDB_CHAIN_TO_UNIPROT_CLASS(PDBS)
print "WITH UNIPROTS", len(TEST.with_uniprots()), TEST.with_uniprots()
print "WITHOUT UNIPROTS", len(TEST.without_uniprots()), TEST.without_uniprots()