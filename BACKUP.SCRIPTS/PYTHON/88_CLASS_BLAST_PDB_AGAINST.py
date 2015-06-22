#BLAST_PDB_AGAINST CLASS #08/02/13
#Depends on properly installed BLAST+ standalone
#Takes PDB|X and blasts it against given uniprot or the entire swissprot if nothing given. Returns uniprot if any that pass cutoff
#filters given. If no uniprot found that met requirements then "None" will be returned.
#WARNING - Make sure PDB|X exists or it will cause an error.
#Coverage is calculated with respect to uniprot, as in, how much does the pdb sequence covers the uniprot sequence.
#A function was implemented to account for mature chains only if any in the uniprot flat files  
 
class BLAST_PDB_AGAINST_CLASS: #Takes PDB|X and blasts it against subject or SWISSPROT database
    
    def __init__ (self, PDBX, IDENTITY_CUTOFF, COVERAGE_CUTOFF, LOG_FOLDER):
        #Setting variables
        self.PDBX=PDBX
        self.PDB=self.PDBX.split("|")[0].upper()
        self.CHAIN=self.PDBX.split("|")[1].upper()
        self.IC=IDENTITY_CUTOFF
        self.CC=COVERAGE_CUTOFF
        self.LOG_FOLDER=LOG_FOLDER
    
    def TO_SUBJECT(self, UNIPROT_LIST): #UNIPROT_ID STRING #If uniprot list was not provided this will produce an error
        from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, UNIPROT_CHAIN_LIMITS
        import subprocess
        
        #Download PDB
        GET_PDB(self.PDB, self.LOG_FOLDER)
        
        #Get PDB_SEQUENCE and save to fasta file if it exists
        PDB_AA, PDB_NUM=PDB_SEQ(self.LOG_FOLDER+"/"+"pdb"+self.PDB.lower()+".ent", self.CHAIN)
        
        if len(PDB_AA)>0: #To account for improperly annotated or wrong chain since it would produce an empty amino acid sequence
            PDB_AA_FASTA=open(self.LOG_FOLDER+"/PDB.fasta","w")
            PDB_AA_FASTA.write(">"+self.PDBX+"\n"+PDB_AA)
            PDB_AA_FASTA.close()
            
            #Get UNIPROT_SEQ for UNIPROT in list and save all to fasta file
            UNIPROT_AA_FASTA=open(self.LOG_FOLDER+"/UNIPROT.fasta", "w")
            for UNIPROT_ID in UNIPROT_LIST:
                UNIPROT_AA=UNIPROT_SEQ(UNIPROT_ID)
                #CHAIN_LIMITS=UNIPROT_CHAIN_LIMITS(UNIPROT_ID) #TO CORRECT FOR MATURE PROTEINS #USE WHEN LOOKING FOR COVERAGE AGAINST
                                                                                                #UNIPROT (1)
                #UNIPROT_AA=UNIPROT_AA[CHAIN_LIMITS[0]: CHAIN_LIMITS[1]] (1)
                UNIPROT_AA_FASTA.write(">"+UNIPROT_ID+"\n"+UNIPROT_AA+"\n")
            UNIPROT_AA_FASTA.close()    
                
            #Do BLAST and only show highest 
            BLAST=subprocess.check_output(["-c", "blastp -query=%s/PDB.fasta -subject=%s/UNIPROT.fasta -num_alignments=1 " 
                                           "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%(self.LOG_FOLDER, self.LOG_FOLDER)], 
                                          shell=True).splitlines()
            
            #Filter by identity>60% and coverage>80 - HERE COVERAGE=ALIGNMENT LENGTH/UNIPROT LENGTH
            if len(BLAST)>0:
                RESULT=BLAST[0].split()
                IDENTITY=float(RESULT[1])
                #COVERAGE=(float(RESULT[3])/float(RESULT[4]))*100 #USE WHEN LOOKING FOR COVERAGE AGAINST UNIPROT (1)
                COVERAGE=float(RESULT[2])
                
                if IDENTITY>float(self.IC) and COVERAGE>float(self.CC):
                    self.UNIPROT=RESULT[0]
                else:
                    self.UNIPROT="None"
            else:
                self.UNIPROT="None"
        else:
            self.UNIPROT="None"
            
        #Return answer
        return self.UNIPROT
    
    def TO_ALL (self): #UNIPROT_ID STRING 
                        #Since uniprot is not expected it aligns to all swissprot database -  THIS FUNCITON DOES TWO ROUNDS OF BLAST
        from FUNCTIONS import GET_PDB, PDB_SEQ, UNIPROT_SEQ, UNIPROT_CHAIN_LIMITS
        import subprocess
    
        #Download PDB
        GET_PDB(self.PDB, self.LOG_FOLDER)
        
        #Get PDB_SEQUENCE and save to fasta file if it exists
        PDB_AA, PDB_NUM=PDB_SEQ(self.LOG_FOLDER+"/"+"pdb"+self.PDB.lower()+".ent", self.CHAIN)
        
        if len(PDB_AA)>0: #To account for improperly annotated or wrong chain since it would produce an empty amino acid sequence        
            PDB_AA_FASTA=open(self.LOG_FOLDER+"/PDB_ALL.fasta","w")
            PDB_AA_FASTA.write(">"+self.PDBX+"\n"+PDB_AA)
            PDB_AA_FASTA.close()
            
            #IMPORTANT - Do BLAST to get HIGHEST CANDIDATE first 
            BLAST2=subprocess.check_output(["-c", "blastp -query=%s/PDB_ALL.fasta -db=swissprot -num_alignments=1 " 
                                           "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%self.LOG_FOLDER], 
                                          shell=True).splitlines()
            
            ##########USE WHEN DOING COVERAGE AGAIN UNIPROT (1)###########                                   
            #Second blast on candidates if any
            """
            if len(BLAST)>0:
                CANDIDATE_UNIPROT=BLAST[0].split()[0].split("|")[3].split(".")[0]
                
                #Get mature sequence of candidate uniprot
                UNIPROT_CANDIDATE_AA=UNIPROT_SEQ(CANDIDATE_UNIPROT)
                CANDIDATE_CHAIN_LIMITS=UNIPROT_CHAIN_LIMITS(CANDIDATE_UNIPROT)
                UNIPROT_CANDIDATE_AA=UNIPROT_CANDIDATE_AA[CANDIDATE_CHAIN_LIMITS[0]: CANDIDATE_CHAIN_LIMITS[1]]
                
                #Save to fasta:
                UNIPROT_CANDIDATE_AA_FASTA=open(self.LOG_FOLDER+"/UNIPROT_CANDIDATE.fasta", "w")
                UNIPROT_CANDIDATE_AA_FASTA.write(">"+CANDIDATE_UNIPROT+"\n"+UNIPROT_CANDIDATE_AA+"\n")
                UNIPROT_CANDIDATE_AA_FASTA.close()
                
                #Blast
                BLAST2=subprocess.check_output(["-c", "blastp -query=%s/PDB_ALL.fasta -subject=%s/UNIPROT_CANDIDATE.fasta -num_alignments=1 "
                                                 "-use_sw_tback=True -outfmt '6 sseqid pident qcovs length slen'"%(self.LOG_FOLDER, self.LOG_FOLDER)],
                                                shell=True).splitlines() 
            
            """
            ################################################################
            
            #Filter by identity>60% and coverage>80 - HERE COVERAGE=ALIGNMENT LENGTH/UNIPROT LENGTH
            if len(BLAST2)>0:
                RESULT=BLAST2[0].split()
                IDENTITY=float(RESULT[1])
                #COVERAGE=(float(RESULT[3])/float(RESULT[4]))*100 #USE WHEN DOING UNIPROT COVERAGE (1)
                COVERAGE=float(RESULT[2])
                print COVERAGE
                
                if IDENTITY>float(self.IC) and COVERAGE>float(self.CC):
                    #self.UNIPROT=RESULT[0] #USE WHEN DOING COVERAGE AGAINST UNIPROT (1)
                    self.UNIPROT=RESULT[0].split("|")[3].split(".")[0]
                else:
                    self.UNIPROT="None"
            else:
                self.UNIPROT="None"
            
            """#USE WHEN DOING COVERAGE AGAINST UNIPROT (1)
            else:
                self.UNIPROT="None"
            """
        else:
            self.UNIPROT="None"
        #Return answer
        return self.UNIPROT      
        
TEST=BLAST_PDB_AGAINST_CLASS("9PAP|a", 60, 80, "LOGS")
A=TEST.TO_SUBJECT([ "Q6GZX4", "P07445","P00784", "Q6GZW9"])
B=TEST.TO_SUBJECT(["Q6GZW9"])
print A
print B

C=TEST.TO_ALL()
print C