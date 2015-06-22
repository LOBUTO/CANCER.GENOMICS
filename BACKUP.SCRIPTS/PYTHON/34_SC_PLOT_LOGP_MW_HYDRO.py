#SPEARMAN AND PEARSON CORRELATION OF LOGP AND MW

import pylab as pl
import scipy.stats
import numpy as np
import subprocess


#LOGP CALCULATIONS
"""
#Including all
LOGP=np.loadtxt("SMILES/32_LOG_P2_SUM") #load data in text file as array
LOGP_SC=scipy.stats.spearmanr(LOGP[:,1],LOGP[:,0])
print LOGP_SC

PLOT1, =pl.plot(LOGP[:,1],LOGP[:,0],"ro")
pl.title("DELTA LIG LOGP VS DELTA SDPA LOGP (ALL)")
pl.xlabel("DELTA SDPA LOGP")
pl.ylabel("DELTA LIG LOGP")
pl.grid(True)
pl.legend([PLOT1], (LOGP_SC),"best")
pl.show()

#Not including zeroes
LOGP2=open("SMILES/32_LOG_P2_SUM")
LOGP2A=[]
for i in LOGP2:
    LOGP2A.append(i.split())
LOGP2A=filter(lambda x: x[0]!="0.0" and x[1]!="0.0",LOGP2A)
LOGP2X=[]
LOGP2Y=[]
for i in LOGP2A:
    LOGP2X.append(i[1])
    LOGP2Y.append(i[0])

LOGP_SC=scipy.stats.spearmanr(LOGP2X,LOGP2Y)
print LOGP_SC

PLOT2, =pl.plot(LOGP2X,LOGP2Y,"bo")
pl.title("DELTA LIG LOGP VS DELTA SDPA LOGP (NON-ZEROES)")
pl.xlabel("DELTA SDPA LOGP")
pl.ylabel("DELTA LIG LOGP")
pl.grid(True)
pl.legend([PLOT2], (LOGP_SC),"best")
pl.show()

#MW CALCULATIONS

#Including zeroes
MW=np.loadtxt("SMILES/32_MW2_SUM")
MW_SC=scipy.stats.spearmanr(MW[:,1],MW[:,0])
print MW_SC

PLOT3, =pl.plot(MW[:,1],MW[:,0],"ro")
pl.title("DELTA LIG MW VS DELTA SDPA MW (ALL)")
pl.xlabel("DELTA SDPA MW")
pl.ylabel("DELTA LIG MW")
pl.grid(True)
pl.legend([PLOT3], (MW_SC),"best")
pl.show()

#Not including zeroes
MW2=open("SMILES/32_MW2_SUM")
MW2A=[]
for i in MW2:
    MW2A.append(i.split())
MW2A=filter(lambda x: x[0]!="0.0" and x[1]!="0.0",MW2A)
MW2X=[]
MW2Y=[]
for i in MW2A:
    MW2X.append(i[1])
    MW2Y.append(i[0])

MW_SC=scipy.stats.spearmanr(MW2X,MW2Y)
print MW_SC

PLOT4, =pl.plot(MW2X,MW2Y,"bo")
pl.title("DELTA LIG MW VS DELTA SDPA MW (NON-ZEROES)")
pl.xlabel("DELTA SDPA MW")
pl.ylabel("DELTA LIG MW")
pl.grid(True)
pl.legend([PLOT4], (MW_SC),"best")
pl.show()

"""
#INDEX
HYDROS=subprocess.check_output("ls", cwd="SMILES/").splitlines()
print HYDROS
HYDROS=filter(lambda x: "35_MEAN" in x, HYDROS)

print HYDROS

for j in HYDROS:
    FILE=open("SMILES/%s"%j)
    FILEA=[]
    for i in FILE:
        FILEA.append(i.split())
    print FILEA
    print len(FILEA)
    
    FILEA_Z=filter(lambda x:x[1]!="nan",FILEA)
    print len(FILEA_Z)
    INDEX_SC=scipy.stats.spearmanr(np.array(FILEA_Z)[:,1],np.array(FILEA_Z)[:,0])
    print INDEX_SC
    PLOT5,=pl.plot(np.array(FILEA_Z)[:,1],np.array(FILEA_Z)[:,0],"g*")
    pl.title("DELTA_LIG_VS_DELTA_SDPA_MEAN_%s(ALL)"%j[8:])
    pl.xlabel("DELTA SDPA")
    pl.ylabel("DELTA LIG")
    pl.grid(True)
    pl.legend([PLOT5],(INDEX_SC),"best")
    pl.savefig("PLOTS/DELTA_LIG_VS_DELTA_SDPA_MEAN_%s(ALL)"%j[8:])
    pl.close()
    
    FILEA_NOZ=filter(lambda x:x[0]!="0.0" and x[1]!="nan" and x[1]!="0.0",FILEA)
    print len(FILEA_NOZ)
    INDEX_SC=scipy.stats.spearmanr(np.array(FILEA_NOZ)[:,1],np.array(FILEA_NOZ)[:,0]) #np.array can be used to make 2D list into array
    print INDEX_SC
    PLOT6,=pl.plot(np.array(FILEA_NOZ)[:,1],np.array(FILEA_NOZ)[:,0],"b*")
    pl.title("DELTA_LIG_VS_DELTA_SDPA_MEAN_%s(NON_ZEROES)"%j[8:])
    pl.xlabel("DELTA SDPA")
    pl.ylabel("DELTA LIG")
    pl.grid(True)
    pl.legend([PLOT6],(INDEX_SC),"best")
    pl.savefig("PLOTS/DELTA_LIG_VS_DELTA_SDPA_MEAN_%s(NON_ZEROES)"%j[8:])
    pl.close()





