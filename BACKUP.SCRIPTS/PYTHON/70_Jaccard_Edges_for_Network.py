#METABOLITE NETWORK WITH PROTEIN JACCARD AS EDGE WEIGHTS - PREPARED WEIGHTED FILE BASED ON JACCARD
def JACCARD(set1, set2): #Feed 2 sets of elements and returns float coefficient 
    result=float(len(set1&set2))/float(len(set1|set2))  
    return result

FILE_IN=open("NETWORK/ENDOGENOUS.ADJ")
END_ADJ=FILE_IN.read().splitlines()
FILE_IN.close()
print END_ADJ[:2]

print len(set(END_ADJ[0].split()[1:]) | set(END_ADJ[1].split()[1:]))
print set(END_ADJ[0].split()[1:]) | set(END_ADJ[1].split()[1:])

print len(set(END_ADJ[0].split()[1:]) & set(END_ADJ[1].split()[1:])) / len(set(END_ADJ[0].split()[1:]) | set(END_ADJ[1].split()[1:]))

#Filter zero degree nodes into file and keep non_zero nodes into list
FILE_OUT1=open("NETWORK/DIRECTED_NETWORK/ZERO_DEGREE_NODES", "w")
NON_ZERO=[]
for line in END_ADJ:
    if len(line.split())==1:
        FILE_OUT1.write(line+"\n")
    else:
        NON_ZERO.append(line)
FILE_OUT1.close()
print len(NON_ZERO)

#Do Jaccard and store weights into files
FILE_OUT2=open("NETWORK/DIRECTED_NETWORK/NODE_EDGE_WEIGHTS","w")
for met1 in NON_ZERO:
    for met2 in NON_ZERO:
        if met1!=met2 and JACCARD(set(met1.split()[1:]), set(met2.split()[1:]))>0:
            PRE_OUT=sorted([met1.split()[0],met2.split()[0]]) +[str(JACCARD(set(met1.split()[1:]), set(met2.split()[1:])))]
            FILE_OUT2.write(PRE_OUT[0]+" "+PRE_OUT[1]+" "+PRE_OUT[2]+"\n")
FILE_OUT2.close()

#Post-processed in terminal for uniqueness, obtained 184212 edges
#DONE