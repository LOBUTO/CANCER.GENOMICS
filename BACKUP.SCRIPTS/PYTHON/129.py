#GENE_LENGTH PROCESSING FROM UNIPROT
#01/20/14

import pickle

#Load file
FILE_IN1=open("DATABASES/UNIPROT/UNIPROT_HOMO_TAB.tab")
UNIPROT=[ [ X.split("\t")[4], X.split("\t")[6]] for X in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()

#Get rid of empty genes
print len(UNIPROT)
UNIPROT=filter(lambda x: x[0]!="", UNIPROT)
print len(UNIPROT)

#Make dictionary - KEEP IN MIND UPPER() WHEN USING IT!!!!
UNIPROT_DICT=dict((X.strip().upper(),int(Y[1])) for Y in UNIPROT for X in Y[0].replace(";","").split())

#012014-UPDATE WITH UNREVIEWD DATA
UNIPROT_DICT['COMMD3-BMI1']=469
UNIPROT_DICT["GCOM1"]=765


#Pickle out
PICKLE_OUT1=open("DATABASES/UNIPROT/OBJECTS/012014_UNIPROT_GENE_LENGTH.pi", "w")
pickle.dump(UNIPROT_DICT, PICKLE_OUT1)
PICKLE_OUT1.close()
