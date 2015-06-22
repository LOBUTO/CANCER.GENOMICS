#Format i.ALN in capra folder
ALN=open("capra_bioinf08_data/1.1.1.1_group0_PF00106_30-85~1.1.1.35_group0_PF00106_30-85.aln")
ALNA=[]
for i in ALN:
    ALNA.append(i)
for i in ALNA:print i #List from file

print 6
ALNB=[]
for i in ALNA:
    ALNB.append(i.split()) #make list of space separated strings

for i in ALNB:print i

print 6
ALNC=[]
for i in ALNB:
    if all(i[0:1]!=j for j in ALNC):
        ALNC.append(i[0:1]) #makes list of unique i[0:1]

for i in ALNC:print i
print len(ALNC)

print 7
ALND=[]
for i in ALNC:
    ALNE=[] #create empty list to store joined string sequences everytime
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
            
            

    