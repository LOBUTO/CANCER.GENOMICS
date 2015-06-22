#PARSE AND CALCULATE DELTA_LOGP and DELTA_MW pairs of PAIR_L_SDPA

import ast #useful for analysis list form in string
import subprocess
import numpy

def unique(LIST):
    DUMMY=[]
    for STRING in LIST:
        if all(STRING!=DE for DE in DUMMY):
            DUMMY.append(STRING)
    LIST=DUMMY 
    return LIST

PAIR_L_SDPA=open("LIGAND_OUTPUT/PAIR_L_SDPA_NO2").readlines() #276 TOTTAL
print PAIR_L_SDPA
print len(PAIR_L_SDPA)

PAIR_L_SDPA=filter(lambda x: x[0]=="[",PAIR_L_SDPA) #Filtering for those that have no sdpA, those lines start with a digit
                                                    #instead of "["
                                                    #159 TOTAL
print PAIR_L_SDPA
print len(PAIR_L_SDPA)

#GET RID OF THOSE THAT HAVE "NO COMMON LIGAND" IN AT LEAST ONE FROM PAIR
PAIR_L_SDPA=filter(lambda y: ast.literal_eval(y)[1]!="No common ligand" #ast.literal_eval() converts list embedded in string
                                                                            #into expected list
                                                                            #here used to remove only
                      and ast.literal_eval(y)[2]!="No common ligand", PAIR_L_SDPA) #39 TOTAL 
print PAIR_L_SDPA
print len(PAIR_L_SDPA)

#REMOVE LABEL
PAIR=[]
for i in PAIR_L_SDPA:
    i=ast.literal_eval(i)
    PAIR.append(i[1:])
print PAIR

#REMOVE METALS AND IONS
PAIR_U=[]
for i in PAIR: #FOR PAIR
    PAIR_L=[]
    for j in i: #FOR EACH GROUP IN PAIR
        j=filter(lambda x: len(x)>2 and x!="SO4", j)
        PAIR_L.append(j)
    PAIR_U.append(PAIR_L)
print len (PAIR_U)

#FILTER THOUSE THAT HAVE LOST ALL LIGANDS DUE TO REMOVAL OF IONS
PAIR_U=filter(lambda z: len(z[0])>1 and len(z[1])>1,PAIR_U) #38 TOTAL
print len(PAIR_U)

######CONSTRUCT#########

#GET ONES THAT HAVE MORE THAN ONE LIGAND PER GROUP - WORK ON THEM
LIGN_PLUS=[]
for l in PAIR_U:
    for k in l:
        if len(filter(lambda w: len(w)<4,k))>1:
            LIGN_PLUS.append(l)

LIGN_PLUS=unique(LIGN_PLUS) #TOTAL 7
print LIGN_PLUS

#REMOVE THOSE THAT HAVE MORE THAN ONE FROM MAIN LIST
for l in LIGN_PLUS:
    PAIR_U.remove(l) #TOTAL 31 NOW

#MAKE ONE LIGAND TO RESIDUES LISTS
PLUS_GROUPS=[]
for l in LIGN_PLUS:
    GROUP1=l[0]
    GROUP2=l[1]
    GROUP1_L=filter(lambda z: len(z)<4, GROUP1) #Ligands in group 1
    GROUP1_R=filter(lambda x: len(x)>4, GROUP1) #Residues in group 1
    GROUP2_L=filter(lambda y: len(y)<4, GROUP2) #Ligands in group 2
    GROUP2_R=filter(lambda w: len(w)>4, GROUP2) #Residues in group 2
    
    #From group1
    GROUPS1=[]
    for m in GROUP1_L:
        GROUPS1.append(m.split() + GROUP1_R) #Add residues to each ligand in group 1

    #From group2
    GROUPS2=[]
    for n in GROUP2_L:
        GROUPS2.append(n.split()+GROUP2_R) #Add residues to each ligand in group 2

    #Combine
    for o in GROUPS1:
        for p in GROUPS2:
            PLUS_GROUPS.append([o,p]) #TOTAL 18

#ADD MAIN GROUP TO PLUS_GROUPS
ALL_GROUPS=PAIR_U+PLUS_GROUPS #NOW 49
print ALL_GROUPS

#IMPORT SMILES
SMILES=open("Components-smiles-oe.smi").readlines()


LOG_P=open("SMILES/32_LOG_P2_SUM", "w")
MW=open("SMILES/32_MW2_SUM","w")
#CALCULATIONS
for r in ALL_GROUPS:
    CALC1=r[0]
    CALC2=r[1]
    LIG1=CALC1[0]
    LIG2=CALC2[0]
    
    #Get SMILES for first ligand
    for s in SMILES:
        if s.split()[1]==LIG1:
            print s
            FILE_SMI1=open("SMILES/32_smi/%s.smi"%"".join(CALC1),"w")
            FILE_SMI1.writelines(" ".join(s.split()))
            FILE_SMI1.close()
    CALC1A=subprocess.Popen(["/usr/local/bin/obprop", "SMILES/32_smi/%s.smi"%"".join(CALC1)],stdout=subprocess.PIPE)
    CALC1B=CALC1A.communicate()[0].splitlines()
    
    for t in CALC1B:
        if t[0:4]=="logP":
            LIG1_LOGP=float(t.split()[1]) #LOGP OF FIRST LIGAND
        elif t[0:10]=="mol_weight":
            LIG1_MW=float(t.split()[1]) #MW OF FIRST LIGAND
    
    #Get SMILES for second ligand
    for ss in SMILES:
        if ss.split()[1]==LIG2:
            print ss
            FILE_SMI2=open("SMILES/32_smi/%s.smi"%"".join(CALC2),"w")
            FILE_SMI2.writelines(" ".join(ss.split()))
            FILE_SMI2.close()
    CALC2A=subprocess.Popen(["/usr/local/bin/obprop", "SMILES/32_smi/%s.smi"%"".join(CALC2)],stdout=subprocess.PIPE)
    CALC2B=CALC2A.communicate()[0].splitlines()
    
    for tt in CALC2B:
        if tt[0:4]=="logP":
            LIG2_LOGP=float(tt.split()[1]) #LOGP OF SECOND LIGAND
        elif tt[0:10]=="mol_weight":    #MW OF SECOND LIGAND
            LIG2_MW=float(tt.split()[1])
    
    #Get SUM SMILES for first residues
    RES1=CALC1[1:]
    RES1_LOGP=[]
    RES1_MW=[]
    for u in RES1:
        RES1_S=u.split("-")[1]
        for sss in SMILES:
            if sss.split()[1]==RES1_S:
                print sss
                FILE_SMI3=open("SMILES/32_smi/%s.smi"%RES1_S, "w")
                FILE_SMI3.writelines(" ".join(sss.split()))
                FILE_SMI3.close()
        CALC1C=subprocess.Popen(["/usr/local/bin/obprop","SMILES/32_smi/%s.smi"%RES1_S],stdout=subprocess.PIPE)
        CALC1D=CALC1C.communicate()[0].splitlines()
        
        for v in CALC1D:
            if v[0:4]=="logP":
                RES1_LOGP.append(float(v.split()[1]))
            elif v[0:10]=="mol_weight":
                RES1_MW.append(float(v.split()[1]))
    
    RES1_LOGP_MEAN=sum(RES1_LOGP)  #LOGP AVERAGE OF FIRST RESIDUE
    RES1_MW_MEAN=sum(RES1_MW)    #MW AVERAGE OF FIRST RESIDUE
    

    #Get SUM SMILES of second residues
    RES2=CALC2[1:]
    RES2_LOGP=[]
    RES2_MW=[]
    for x in RES2:
        RES2_S=x.split("-")[1]
        for ssss in SMILES:
            if ssss.split()[1]==RES2_S:
                print ssss
                FILE_SMI4=open("SMILES/32_smi/%s.smi"%RES2_S, "w")
                FILE_SMI4.writelines(" ".join(ssss.split()))
                FILE_SMI4.close()
        CALC2C=subprocess.Popen(["/usr/local/bin/obprop","SMILES/32_smi/%s.smi"%RES2_S],stdout=subprocess.PIPE)
        CALC2D=CALC2C.communicate()[0].splitlines()
        
        for y in CALC2D:
            if y[0:4]=="logP":
                RES2_LOGP.append(float(y.split()[1]))
            elif y[0:10]=="mol_weight":
                RES2_MW.append(float(y.split()[1]))
    
    RES2_LOGP_MEAN=sum(RES2_LOGP)  #LOGP AVERAGE OF SECOND RESIDUE
    RES2_MW_MEAN=sum(RES2_MW)    #MW AVERAGE OF SECOND RESIDUE
    
    #Get DELTAs (FIRST-SECOND)
    DELTA_LIG_LOGP=LIG1_LOGP-LIG2_LOGP
    DELTA_RES_LOGP=RES1_LOGP_MEAN-RES2_LOGP_MEAN
    DELTA_LIG_MW=LIG1_MW-LIG2_MW
    DELTA_RES_MW=RES1_MW_MEAN-RES2_MW_MEAN
    
    LOG_P.write(str(DELTA_LIG_LOGP)+" "+str(DELTA_RES_LOGP)+"\n")
    MW.write(str(DELTA_LIG_MW)+" "+str(DELTA_RES_MW)+"\n")
    
LOG_P.close()
MW.close()