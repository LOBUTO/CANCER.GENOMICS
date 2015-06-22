#5-MODIFIED
#i.aln, i.ligdist5 AND i.sdpo lists from folders + X-over to get M1_O 

import glob

ALN=[]
for name in glob.glob("capra_bioinf08_data/*.aln"): #glob.glob (*.*) returns a list of pathnames that match pathname specified
    print name
    ALN.append(name)

print ALN
for i in ALN: print i

LIGD5=[]
for name in glob.glob("capra_bioinf08_data/*.ligdist5"): #glob.glob (*.*) returns a list of pathnames that match pathname specified
    print name
    LIGD5.append(name)

print LIGD5
for i in LIGD5: print i

ALNR=[]
for i in ALN:
    ALNR.append(i[:-4])#to reduce the name without extension to X-ref to LIGD5R

print ALNR

print 6
LIGD5R=[]
for i in LIGD5:
    LIGD5R.append(i[:-9])#to reduce the name without extension to X-ref to ALNR
    
print LIGD5R

#MODIFICATION
SDPO=[]
for name in glob.glob("capra_bioinf08_data/*.sdpO"): #storing the i.sdpO lists to add later to M1
    print name
    SDPO.append(name)
    
SDPOR=[]
for i in SDPO:
    SDPOR.append((i[:-5])) #reducing the extension to match
    

M1=[] #FIRST MATRIX TO STORE ['ID_X','Pfam','EC','Sequence','W'] ('W' is the first residue within 5A of ligand with respect to i.aln starting with 0)
for i in ALNR:
    for j in LIGD5R:
        for k in SDPOR:
            if i==j==k:#if the names are the same THEN THEY MUST COME FROM THE SAME SET, so during each call it will append to M1 only from the same set!
                
                #MAKE SURE TO OPEN BOTH FROM THE BEGINNING OTHERWISE IF "j" or "k" ARE CALLED LATER IT WILL CONFUSE WITH ANOTHER CALLING
                i=i+".aln" #add proper extension to i.aln files
                ALNG=open(i) #For list i
                j=j+".ligdist5" #add proper extension to j.ligdist5 files
                LD5C=open(j) #For list j
                k=k+".sdpO" #add proper extension to k.sdpO files
                SDPOA=open(k) #For list k
                
                ###############
                #Format each i.ALN in capra folder that is equal to i.ligdist5 and i.spdO in capra folder
                ALNA=[] #RESTART
                for i in ALNG:
                    ALNA.append(i)
                for i in ALNA:print i #List from file
               
                ALNB=[]
                for i in ALNA:
                    ALNB.append(i.split()) #make list of space separated strings
                
                for i in ALNB:print i
                
                ALNC=[]
                for i in ALNB:
                    if all(i[0:1]!=j for j in ALNC):
                        ALNC.append(i[0:1]) #makes list of unique i[0:1]
                
                for i in ALNC:print i
                print len(ALNC)
                
                ALND=[]
                for i in ALNC:
                    ALNE=[] #create empty list to store joined string sequences every time
                    for j in ALNB:        
                        if i==j[0:1]:
                            ALNE=ALNE+j[1:2] #if they have the same tag then join all the sequence strings into a list of strings called ALNE
                    ALNE="".join(ALNE).split() #join list of sequence strings into a single string, however we need to make it into a list with .split() to use append later         
                    ALND.append(i+ALNE) #spit unique tag + whole sequence
                
                for i in ALND:print i
                print len (ALND)
                
                print 8
                ALNF=[] #need to make one more formating to separate by ID_X so I can call from within i.ligdist5 later on
                for i in ALND:
                    for j in i[0:1]: #because i[0:1] is not a string but a list element, so we need to call the string within it first
                        ALNF.append(j.split("|")+i[1:2])
                
                for i in ALNF: print i 
                print len(ALNF)
                
                #Format i.ligdist5 from capra folder
                LD5A=[] #RESTART
                for i in LD5C:
                    LD5A.append(i.split()) #i.split() produces a list of the things being splitted, append in this case is appending lists
                
                print LD5A
                
                LD5B=[]
                for i in LD5A:
                    for j in i[0:1]:
                        if (j!="#" and j!="*"):
                            LD5B.append(i) #get only ID_X + all residues within 5A of ligand (with respect to i.aln starting at 0)
                        
                
                for i in LD5B: print i
                
                #Format i.sdpO from capra folder
                SDPOB=SDPOA.readlines() #separate into strings
                SDPOC=SDPOB[0].split() #store and split first line that contains name of file set
                SDPOD=SDPOC[1:2] #Keep name of file set into list format
                    
                #X-REF BOTH LISTS WHICH COME FROM THE SAME SET BY ID_X + SDPOC OF SAME FILE SET
                for i in LD5B:
                    for j in ALNF:
                        if i[0:1]==j[0:1]: #if any ID_X is the same in both i.aln and j.ligdist5, then:
                            
                            #CONSTRUCT M1 MATRIX
                            M1.append(j[0:1]+j[1:2]+j[2:3]+j[3:4]+i[1:2]+SDPOD) #using first number of ALNF

M1.sort(key=lambda X:X[0:1]) #Sorts M1 based on first element in each list            
for i in M1: print i
print len(M1) #643

#Look for redundancies
M1U=[]
for i in M1:
    if all(i!=j for j in M1U):        
        M1U.append(i) #M1 SORTED AND UNIQUE
    else:
        print i
print len(M1U) #643

output=open("M1_O","w") #MODIFIED TO M1_O
for i in M1U:
    output.write("%s\n"%" ".join(i))


    



