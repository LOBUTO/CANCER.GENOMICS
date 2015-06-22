#Getting unique files out of pdb_seqres_AUNIQUE and M1_NO1THG to CLUSTALO/M1 and CLUSTALO/PDAU

M1T=open("M1_O_NO1THG") #Open up M1_O_NO1THG
M1TA=M1T.readlines()
for i in M1TA:
    if i[5]=="_":
        M1TC=i.split()
        PRE_ID_X=M1TC[0]
        ID_X=PRE_ID_X[0:5]+"A"
        SET=M1TC[5]
    
    else:
        M1TC=i.split()
        ID_X=M1TC[0]
        SET=M1TC[5]
    
    FILE=ID_X+"_"+SET
    M1TB=open("CLUSTALO/M1/%s"%FILE,"w")  #641, it matches
    M1TB.writelines(i) #Write file with elements of each single line
    M1TB.close()
        
        
SRU=open("SEQRES/pdb_seqres_AUNIQUE") #Open up AUNIQUE
SRUA=SRU.readlines()
for j in SRUA:
    SRUB=open("CLUSTALO/PDAU/%s"%j[0:6],"w") #235, it matches
    SRUB.writelines(j)
    SRUB.close()
    
    

