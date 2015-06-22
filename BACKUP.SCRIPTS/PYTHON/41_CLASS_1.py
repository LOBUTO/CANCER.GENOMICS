# SCRIPT 1 - BUILDING GO_ID_PFAM DICTIONARIES AND GO CLUSTERS

def unique(LIST): #LIST - FILTER FOR UNIQUENESS IN LIST
    DUMMY=[]
    for STRING in LIST:
        if STRING not in DUMMY:
            DUMMY.append(STRING)
    return DUMMY

def MIX(LIST1, LIST2,FIX=0): #LIST - MIXES EVERY CONTENT OF LIST1(n) WITH EVERY ELEMENT OF LIST2(m) RETURNING LIST OF mxn ELEMENTS
                            #IF FIX=1 THEN JOINS EVERY ELEMENT OF LIST 2 WITH ALL ELEMENTS OF LIST 1 TO MAKE nxLIST1 LISTS
                            #IF FIX=2 THEN JOINS EVERY ELEMENT OF LIST 1 WITH ALL ELEMENTS OF LIST 2 TO MAKE LIST2xm LISTS
    LIST=[]
    if FIX==0:
        for i in LIST1:
            for j in LIST2:
                LIST.append([i,j])
    elif FIX==1:
        for l in LIST2:
            LIST.append(LIST1+[l])
    elif FIX==2:
        for k in LIST1:
            LIST.append([k]+LIST2)
    
    return LIST

"""FILE_IN=open("CLASS_PROJECT/gene_ontology_ext.obo").readlines() #Parsing gene_ontology classification file


FILE_IN=FILE_IN[25:-49]
print FILE_IN[0:10]
print FILE_IN[8]
print len(FILE_IN[8])


FILE_OUT=open("CLASS_PROJECT/GO_PRE_DICT","w")

X=0
while X!=len(FILE_IN):
    if FILE_IN[X][0:2]=="id":
        FILE_OUT.write(FILE_IN[X].split()[1]+" ")
        X=X+1
    elif FILE_IN[X][0:9]=="namespace":
        FILE_OUT.write(FILE_IN[X].split()[1]+" ")
        X=X+1
    elif FILE_IN[X][0:4]=="is_a": 
        FILE_OUT.write(FILE_IN[X].split()[1]+" ")
        X=X+1
    
    elif FILE_IN[X]=="\n":
        FILE_OUT.write("\n")
        X=X+1
    
    else:
        X=X+1

FILE_OUT.close()
"""

def GO_DICT_LIST(GO_PRE_FILE, OUT_FILE_DIC): #FILE - First file should be in the format "GO, gene ontology, GO(from is_a)
                                        #Outputs FILE of lines GO(is_a) " " GO
                                        #WARNING - MAKE SURE TO MAKE UNIQ
    FILE=open(GO_PRE_FILE)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split())    
    
    #Get rid of obsolete GOs 
    FILEA=filter(lambda x: len(x)>2,FILEA)
        
    #Get is_a TAGS:
    TAGS=[]
    for j in FILEA:
        TAGS.extend(j[2:])
    
    #MAKE DICT:
    FILEB=open(OUT_FILE_DIC,"w")
    for h in TAGS:
        for g in FILEA:
            if h in g[2:]:
                FILEB.write(h+" "+g[0]+"\n")

    
    FILEB.close()

#TEST SPEED OF DICT

def MAKE_DICT_LIST(FILE_IN): #DICTIONARY - Makes GO Term dictionary with OUT_FILE_DIC
    FILE=open(FILE_IN)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split())
    
    print FILEA
    
    MOL_DICT={} #VERY FAST - DICTIONARY
    for key,val in FILEA:
        MOL_DICT.setdefault(key,[]).append(val)
    
    FILE.close()
    return MOL_DICT
    
def END_GO(GO_DICT): #Takes GO dictionary from MAKE_DICT_LIST and returns a list of terms that have no children
    DICT_VALUES=[]
    for i in GO_DICT.values():
        for j in i:
            DICT_VALUES.append(j)

    END_GO=[]
    for i in DICT_VALUES:
        if GO_DICT.has_key(i)==False:
            END_GO.append(i)

    return END_GO
    
SEPARATOR=[]

#GET GO TERM DICTIONARY
MOL_DICT=MAKE_DICT_LIST("CLASS_PROJECT/DICT/GO_DICT_MOL")

#CHECK FOR GO TERMS THAT HAVE TWO TERMINAL GO TERMS
#First, get all terms that have no children
END_GO=END_GO(MOL_DICT)

#Second get parents of terminal GO terms
END_KEYS=[]
for key, val in MOL_DICT.iteritems():
    for i in END_GO:
        if i in val:
            END_KEYS.append(key)

END_KEYS=unique(END_KEYS) #1688

#Third, filter for end terms that only have two children GO terms
END_KEYS= filter(lambda x: len(MOL_DICT[x])==2 , END_KEYS) #391 GO TERMS

#GET GO_PFAM DICTIONARY
GO_PFAM_DICT=MAKE_DICT_LIST("CLASS_PROJECT/GO_PFAM")

#DEVELOP GROUPS AND SUBGROUPS 

#First, get all children GO Terms # GOT 755 unique(parsed in terminal)
GO_FOR_DICT=open("CLASS_PROJECT/GO_FOR_DICT", "w")
for i in END_KEYS:
    for j in MOL_DICT[i]:
        GO_FOR_DICT.write(j +"\n")
GO_FOR_DICT.close()

#Second parse in Terminal
    #Got ID_PFAM_FOR_DICT (PDB or SWISS to PFAM) of interest #27033 found    
    #Got GO_ID_FOR_DICT (GO to PDB or SWISS) of interest #27112 unique sorted

#Third build dictionaries GO_IDs and ID_PFAM dictionary
GO_ID=MAKE_DICT_LIST("CLASS_PROJECT/GO_ID_FOR_DICT")
ID_PFAM=MAKE_DICT_LIST("CLASS_PROJECT/ID_PFAM_FOR_DICT")

#Third get subgroup families #GOT 218 GROUPS
X=1
for i in END_KEYS:
    #Get GO terms of both children
    GO1=MOL_DICT[i][0]
    GO2=MOL_DICT[i][1]
    
    if GO_ID.has_key(GO1) and GO_ID.has_key(GO2):
        #Get the IDs for each subgroup
        GO1_ID=GO_ID[GO1]
        GO2_ID=GO_ID[GO2]
        
        #Get the PFAMs from each ID in each subgroup
        GO1_PFAM=[]
        for ID1 in GO1_ID:
            if ID_PFAM.has_key(ID1):
                GO1_PFAM.extend(ID_PFAM[ID1])
        GO2_PFAM=[]
        for ID2 in GO2_ID:
            if ID_PFAM.has_key(ID2):
                GO2_PFAM.extend(ID_PFAM[ID2])
        print GO1_PFAM
        print GO2_PFAM
            
        #Compare them and save equals as groups of subgroups
        for PFAM1,PFAM2 in MIX(GO1_PFAM,GO2_PFAM):
            if PFAM1==PFAM2:
                FILE_OUT=open("CLASS_PROJECT/GROUPS/%s"%(str(X)+"_"+GO1[3:]+"_"+PFAM1+"~"+GO2[3:]+"_"+PFAM2),"w")
                for j in GO1_ID:
                    if ID_PFAM.has_key(j) and PFAM1 in ID_PFAM[j]:
                        FILE_OUT.write(GO1+" "+j+" "+PFAM1+" "+"\n")
                for k in GO2_ID:
                    if ID_PFAM.has_key(k) and PFAM1 in ID_PFAM[k]:
                        FILE_OUT.write(GO2+" "+k+" "+PFAM2+" "+"\n")
                FILE_OUT.close()
    
    print X
    X=X+1