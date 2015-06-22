#FIXES FOR LACK OF HET_ID IN METABOLITES_ENDOGENOUS_FIXED

FILE1=open("NETWORK/METABOLITES_ENDOGENOUS").read().splitlines()
FILE2=open("NETWORK/METABOLITES_ENDOGENOUS_FIXED").read().splitlines()

DIFF=[]
for i in FILE1:
    if all(i.split("$")[0].split("|")[0]!=j.split("$")[0].split("|")[0] for j in FILE2):
        DIFF.append(i)
print len(DIFF)
for i in DIFF: print i
