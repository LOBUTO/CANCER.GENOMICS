#FORMAT pdbtosp.text FOR UNIPROT NON-COMPLEX CROSS-REFERENCE #NO COMPLEXES, DEPRECATED
import pickle
from FUNCTIONS import unique

FILE_IN1=open("DATABASES/pdbtosp_201308.txt").read().splitlines()
print len(FILE_IN1)

for record in FILE_IN1:
    if record.split()[0]=="1A02":
        print record.split()
    elif record.split()[0]=="NFAC2_HUMAN":
        print record.split()
    elif record.split()[0]=="1A03":
        print record.split()

#Filter out complexes
FILE_IN1=filter(lambda x: len(x.split()[0])==4 and len(x.split())<7, FILE_IN1)
print len(FILE_IN1)
print 7
for record in FILE_IN1:
    if len(record.split()[-1])!=8: #CHECKS
        print record

print 8
#Filter for human only -OPTIONAL
#FILE_IN1=filter(lambda x: x.split()[-2].split("_")[1]=="HUMAN", FILE_IN1)
print len(FILE_IN1)

#Get PDB_ID and UNIPROT only
FILE_IN1=[[x.split()[0],x.split()[-1].strip("(").strip(")")] for x in FILE_IN1]
print FILE_IN1[:3]

#Convert to Dictionary
DICT={}
for record in FILE_IN1:
    DICT[record[0]]=record[1]
print DICT["10GS"]

#Save as dictionary object #THE "ALL" RECORD WAS REPLACED IN SCRIPT 72.PY WHICH INCLUDES COMPLEXES 
PICKLE_OUT1=open("DATABASES/OBJECTS/DICT_PDB_TO_UNIPROT_ALL.pi", "w")
pickle.dump(DICT, PICKLE_OUT1)
