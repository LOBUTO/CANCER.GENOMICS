#CALCULATE REGRESSION COEFFICIENT - DIFFERENT HYDROPHOBICITY INDEX

import ast
import subprocess
import numpy

def unique(LIST):
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

PAIR=open("LIGAND_OUTPUT/PAIR_L_SDPA_NO2").readlines()
PAIR=filter(lambda x:x[0]=="[",PAIR)
PAIRA=[]
for i in PAIR: 
    PAIRA.append(ast.literal_eval(i)[1:])
print len(PAIRA)

PAIRA=filter(lambda y:y[0]!="No common ligand" and y[1]!="No common ligand", PAIRA)
print len(PAIRA)

print PAIRA

def FILT_LIG(LIST): #LIST - FILTERS FOR LIGANDS IN LIST
    LIST=filter(lambda x:len(x)<4,LIST)
    return LIST
def FILT_RES(LIST): #LIST- FILTER FOR ANYTHING BESIDES LIGANDS IN LIST
    LIST=filter(lambda x:len(x)>4,LIST)
    return LIST    

def MIX(LIST1, LIST2,FIX=0): #MIXES EVERY CONTENT OF LIST1(n) WITH EVERY ELEMENT OF LIST2(m) RETURNING LIST OF mxn ELEMENTS
                            #IF FIX=1 THEN JOINS EVERY ELEMENT OF LIST 2 WITH ALL ELEMENTS OF LIST 1 TO MAKE nxLIST1 LISTS
                            #IF FIX=2 THEN JOINS EVERY ELEMENT OF LIST 1 WITH ALL ELEMENTS OF LIST 2 TO MAKE LIST2xm LISTS
    LIST=[]
    if FIX==0:
        for i in LIST1:
            for j in LIST2:
                LIST.append([i,j])
    elif FIX==1:
        for l in LIST2:
            LIST.append(LIST1+[l])
    elif FIX==2:
        for k in LIST1:
            LIST.append([k]+LIST2)
    
    return LIST
 
def SMI_FILE(MOL,FILE_IN): #SAVES SMILE FILE, MOL IS MOLECULE OF INTEREST AND FILE_IN IS FILE TO STORE IT IN
    SMILES=open("Components-smiles-oe.smi").readlines() #UPDATE THIS IF SMILES DATABASE IS PLACED SOMEWHERE ELSE
    for i in SMILES:
        if i.split()[1]==MOL:
            FILE=open(FILE_IN,"w")
            FILE.writelines(" ".join(i.split()))
            FILE.close()

def LOGP(FILE_IN):
    LOG_P=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    LOG_PA=LOG_P.communicate()[0].splitlines()
    for i in LOG_PA:
        if i[0:4]=="logP":
            LOG_PB=float(i.split()[1])
    return LOG_PB

def MW(FILE_IN):
    MWA=subprocess.Popen(["/usr/local/bin/obprop",FILE_IN],stdout=subprocess.PIPE)
    MWB=MWA.communicate()[0].splitlines()
    for i in MWB:
        if i[0:10]=="mol_weight":
            MWB=float(i.split()[1])
    return MWB

PAIRB=[]
for i in PAIRA:
    FIRST=i[0]
    SECOND=i[1]
    LIG1=FILT_LIG(FIRST)
    RES1=FILT_RES(FIRST)
    LIG2=FILT_LIG(SECOND)
    RES2=FILT_RES(SECOND)
    
    PAIRB.extend(MIX(MIX(LIG1,RES1,FIX=2), MIX(LIG2,RES2,FIX=2)))

print PAIRB

HYDROS=subprocess.check_output("ls", cwd="DATABASES/").splitlines()
if "ALL_HYDRO_SCALES" in HYDROS: HYDROS.remove("ALL_HYDRO_SCALES")
if "HYDRO_INDEX.txt" in HYDROS: HYDROS.remove("HYDRO_INDEX.txt")

print HYDROS 

for o in HYDROS:
    HYDROS_IN=open("DATABASES/%s"%o)
    HYDROA=[]
    for p in HYDROS_IN: HYDROA.append(p.split())

    LOG_PZ=open("SMILES/35_MEAN_%s"%o, "w")
    for i in PAIRB:
        print i
        LIG1=FILT_LIG(i[0][0])
        RES1=FILT_RES(i[0][1:])
        LIG2=FILT_LIG(i[1][0])
        RES2=FILT_RES(i[1][1:])
        
        #LIG1 SMILES
        SMI_FILE(LIG1,"SMILES/35_smi/%s.smi"%LIG1)
        LIG1_LOGP=LOGP("SMILES/35_smi/%s.smi"%LIG1)
        LIG1_MW=MW("SMILES/35_smi/%s.smi"%LIG1)
        
        #LIG2 SMILES
        SMI_FILE(LIG2,"SMILES/35_smi/%s.smi"%LIG2)
        LIG2_LOGP=LOGP("SMILES/35_smi/%s.smi"%LIG2)
        LIG2_MW=MW("SMILES/35_smi/%s.smi"%LIG2)
    #    
        XRES1=[]
        for j in RES1:
            for l in HYDROA:
                if j.split("-")[1]==l[0]:
                    XRES1.append(float(l[1]))
        XRES1M=numpy.mean(XRES1)
        print XRES1M
    #
        
        XRES2=[]
        for k in RES2:
            for m in HYDROA:
                if k.split("-")[1]==m[0]:
                    XRES2.append(float(m[1]))
        XRES2M=numpy.mean(XRES2)
        print XRES2M
        
        DELTA_LIG_LOGP=LIG1_LOGP-LIG2_LOGP
        DELTA_RES_LOGP=XRES1M-XRES2M
        
        LOG_PZ.write(str(DELTA_LIG_LOGP)+" "+str(DELTA_RES_LOGP)+"\n")
    
    LOG_PZ.close()
