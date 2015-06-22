#HMDB METABOLITES PROCESSING 2.1 - HELPER
#08/20/13
#Process FDA DRUG DATABASE to get DRUG NAMES.
#DRUGBANK from HMDB has too many endogenous metabolites listed as drugs in their files

#Load file
FILE_IN=open("DATABASES/FDA/drugsatfda/Product.txt")
DATABASE=FILE_IN.read().splitlines()
FILE_IN.close()

#Get names
DRUG_NAMES=list(set([X.split("\t")[7] for X in DATABASE[1:]]))

#Write to output
FILE_OUT=open("DATABASES/FDA/FDA_NAMES", "w")
for record in DRUG_NAMES:
    FILE_OUT.write(record+"\n")
FILE_OUT.close()

#DELETED LAST EMPTY LINE JUST BECAUSE