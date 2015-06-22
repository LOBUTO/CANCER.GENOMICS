#Parsing of pdb_seqres(from PDBANK) into pdb_seqres_AUNIQUE (Aligned pdb seqres and unique ID_X)

SEQRES=open("SEQRES/pdb_seqres.txt")
SEQRESA=[]
for i in SEQRES:
    SEQRESA.append(i.split("\n")) #readlines would have given me "\n" which would be difficult to manage when joining i and i+1
                                #later, splitting at \n will let me have always the second element in each list as the
                                #sequence

print len (SEQRESA)

SEQRESB=[]

X=0
while X<len(SEQRESA): #To do it until we hit all of the lines in SEQRES 
    SEQRESC=SEQRESA[X]+SEQRESA[X+1] #join the tag at X and its sequence in the following line
#    print SEQRESC #To see what is it doing

    SEQRESB.append(SEQRESC)
    X=X+2 #to increase for every ">" character every other line to account for the pdb name first and the sequence assigned
        #to it

    
#CALL SAVE TO SEQRESD FILE IF NECESSARY
#SEQRESD=open("SEQRES/pdb_seqres_aligned.txt", "w")
#for i in SEQRESB:
#    SEQRESD.write("%s\n"%" ".join(i))

#SEQRESD.close() #Get into the habit of closing files


#NEED TO DO IN SAME SCRIPT FILTER TO LIST OBTAINED FROM ID_X FROM M1_NO1THG
M1N=open("M1_NO1THG")
M1NA=[]
for i in M1N:
    i=i.split()
    M1NA.append(i[0]) #cause only need ID_X (first element of list)
    
M1NB=[]
for j in M1NA:
    for k in SEQRESB:
        if j[5]=="_":#For those ID_X in M1_NO1THG that have an empty chain we add "A", standard for PDBANK files
            ID=j[0:5]+"A" 
        else:
            ID=j
        
        ID_X_PRE=k[0]
        ID_X=ID_X_PRE[1:7] #We call the ID_X from SEQRESB and get only the ID_X without the ">"
        if ID==ID_X: #X-Ref them and append only those that match to the M1_NO1THG list (the ones we need) 
            M1NB.append(ID_X.split()+k[2].split()) #Because we only need the ID_X and the sequence that it accompanies, needs to
                                                    #make them into lists with i.split() to be able to save them later

M1NC=[]
for l in M1NB: #filter for uniqueness as list
    if all (l!=m for m in M1NC):
        print l
        M1NC.append(l)

print len(M1NC) #theoretically should be 234 unique ID_X (Not 236  like in 7...py since we are using a database (M1_NO1THG)
                #that has two less copies of 1THG 

OUTPUT=open("SEQRES/pdb_seqres_AUNIQUE","w")
for n in M1NC:
    OUTPUT.write("%s\n"%" ".join(n))

OUTPUT.close()





    
            
        
    
    
    
    
        

