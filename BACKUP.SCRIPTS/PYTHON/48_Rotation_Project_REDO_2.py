#Get LOGP and MW of SDPA pairs
import numpy as np
import subprocess

def SMI_FILE(MOL,FILE_IN): #SAVES SMILE FILE, MOL IS MOLECULE OF INTEREST AND FILE_IN IS FILE TO STORE IT IN
    SMILES=open("Components-smiles-oe.smi").readlines() #UPDATE THIS IF SMILES DATABASE IS PLACED SOMEWHERE ELSE
    for i in SMILES:
        if i.split()[1]==MOL:
            FILE=open(FILE_IN,"w")
            FILE.writelines(" ".join(i.split()))
            FILE.close()

def LOGP(FILE_IN): #FLOAT - GETS LOGP OF SINGLE SMILES CODE IN GIVEN FILE
    LOG_P=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    LOG_PA=LOG_P.communicate()[0].splitlines()
    LOG_PB=""
    for i in LOG_PA:
        if i[0:4]=="logP":
            LOG_PB=float(i.split()[1])
    return float(LOG_PB)

def MW(FILE_IN): #FLOAT - GETS MW OF SINGLE SMILES CODE IN GIVEN FILE
    MWA=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    MWB=MWA.communicate()[0].splitlines()
    MWC=""
    for i in MWB:
        if i[0:10]=="mol_weight":
            MWC=float(i.split()[1])
    return float(MWC)

SEPARATOR=[]

#Parse File
FILE_IN=open("LIGAND_OUTPUT/LIG_AA_PAIR_ALT_SINGLE_SET").read().splitlines()
print FILE_IN

FILE_OUT_LOGP=open("LIGAND_OUTPUT/LOGP_SINGLE_SET", "w")
FILE_OUT_MW=open("LIGAND_OUTPUT/MW_SINGLE_SET","w")
for i in FILE_IN:
    i=i.split("|")
    FIRST_LIG=i[1].split("_")[0]
    FIRST_AA=i[1].split("_")[1:]
    SECOND_LIG=i[2].split("_")[0]
    SECOND_AA=i[2].split("_")[1:]
    
    #GET LIG DELTAS
    FIRST_SMI_LIG=SMI_FILE(FIRST_LIG, "LOGS/%s.smi"%FIRST_LIG)
    SECOND_SMI_LIG=SMI_FILE(SECOND_LIG, "LOGS/%s.smi"%SECOND_LIG)
    DELTA_LIG_LOGP=LOGP("LOGS/%s.smi"%FIRST_LIG)-LOGP("LOGS/%s.smi"%SECOND_LIG)
    DELTA_LIG_MW=MW("LOGS/%s.smi"%FIRST_LIG)-MW("LOGS/%s.smi"%SECOND_LIG)
    
    #GET AA AVG_DELTAS
    FIRST_AA_LOGP,FIRST_AA_MW=[],[]
    for j in FIRST_AA:
        FIRST_SMI_AA=SMI_FILE(j, "LOGS/%s.smi"%j)
        FIRST_AA_LOGP.append(LOGP("LOGS/%s.smi"%j))
        FIRST_AA_MW.append(MW(("LOGS/%s.smi"%j)))

    SECOND_AA_LOGP,SECOND_AA_MW=[],[]    
    for l in SECOND_AA:
        SECOND_SMI_AA=SMI_FILE(l, "LOGS/%s.smi"%l)
        SECOND_AA_LOGP.append(LOGP("LOGS/%s.smi"%l))
        SECOND_AA_MW.append(MW(("LOGS/%s.smi"%l)))        
    
    DELTA_AA_LOGP=np.mean(FIRST_AA_LOGP)-np.mean(SECOND_AA_LOGP)
    DELTA_AA_MW=np.mean(FIRST_AA_MW)-np.mean(SECOND_AA_MW)
    
    FILE_OUT_LOGP.write(str(DELTA_LIG_LOGP)+" "+str(DELTA_AA_LOGP)+"\n")
    FILE_OUT_MW.write(str(DELTA_LIG_MW)+ " "+str(DELTA_AA_MW)+"\n")

FILE_OUT_LOGP.close()
FILE_OUT_MW.close()