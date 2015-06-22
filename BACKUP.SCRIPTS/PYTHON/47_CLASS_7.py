#SCRIPT 7 - ANALYSIS

import pylab as pl
import math
import numpy as np
import scipy.stats

#FIRST - Get percentage of TRUE POSITIVE/PREDICTED - #SHOWED THAT THERE ARE VERY FEW OF PREDICTED SDPA THAT ARE TRUE SDPA
SDPA_TEST_IN=open("CLASS_PROJECT/RESULTS/SDPA_TEST")
SDPA=[]
for i in SDPA_TEST_IN: SDPA.append(i.split())
SDPA_TEST_IN.close()

X=0
Y=0
TPOS_PRE=[]
for i in SDPA:
    TOTAL=len(i[1:])
    YES=" ".join(i[1:]).count("Y")
    NO=" ".join(i[1:]).count("N")
    X=X+YES
    Y=Y+TOTAL
    TPOS_PRE.append((float(YES)/float(TOTAL))*100)
print X, Y

pl.hist(TPOS_PRE,bins=20)
pl.title("True SDPAs Over All Predictions")
pl.xlabel("Percent Accuracy")
pl.ylabel("Number of Predictions")
pl.grid(True)
pl.show() 

#SECOND - Get True Positives (Retrieved) vs Total positives(Relevant) #GOT A NEGATIVE PEARSON CORRELATION COEFF. SHOWS THAT THERE ARE VERY
                                                                        #FEW CASES IN WHICH ALL TRUE POSITIVES FOUND ARE ALL TOTAL POSITIVES
GOLD_TEST_IN=open("CLASS_PROJECT/RESULTS/GOLD_SD")
GOLD=[]
for i in GOLD_TEST_IN:GOLD.append(i.split())
GOLD_TEST_IN.close()

TOT_TRUE=[]
X=0
while X!=237:
    if len(GOLD[X][1:])>0:
        TOT_POS=len(GOLD[X][1:])
        TRUE_POS=" ".join(SDPA[X][1:]).count("Y")
        TOT_TRUE.append([TOT_POS, (float(TRUE_POS)/float(TOT_POS))])
    X=X+1

PLOT2,=pl.plot(np.array(TOT_TRUE)[:,1], np.array(TOT_TRUE)[:,0], "b*")
INDEX_PC=scipy.stats.spearmanr(np.array(TOT_TRUE)[:,1], np.array(TOT_TRUE)[:,0])
print INDEX_PC
pl.title("Total Positives vs Retrieved Positives")
pl.xlabel("Retrieved Positives")
pl.ylabel("Total Positives")
pl.grid(True)
pl.legend([PLOT2],(INDEX_PC), "best")
pl.show()

#THIRD - Get Precision recall
PRE_REC=[]
X=0
while X!=237:
    if len(GOLD[X][1:])>0:
        TP=" ".join(SDPA[X][1:]).count("Y")
        FP=" ".join(SDPA[X][1:]).count("N")
        
        #Get negatives established by method
        ALL_SDPA=[]
        for i in SDPA[X][1:]:ALL_SDPA.append(int(i.split("-")[0]))    
        
        PFAM=SDPA[X][0].split("|")[1]
        T_N=[]
        for i in range(0,int(PFAM.split("-")[1])-int(PFAM.split("-")[0])):
            if i not in ALL_SDPA:
                T_N.append(i)
        
        TN=0    
        for i in T_N:
            if all(str(i)!=j.split("-")[0] for j in GOLD[X][1:]):
                TN=TN+1
        
        FN=0
        for i in T_N:
            if any(str(i)==j.split("-")[0] for j in GOLD[X][1:]):
                FN=FN+1
        
        RECALL=float(TP)/float(TP+FN)
        PRECISSION=float(TP)/float(TP+FP)
        PRE_REC.append([PRECISSION,RECALL])
    
    X=X+1

PRE_REC=sorted(PRE_REC, key=lambda x : x[1])

PLOT3,=pl.plot(np.array(PRE_REC)[:,1],np.array(PRE_REC)[:,0], 'r-')
pl.title("Precision-Recall Curve - GO-PFAM Data set")
pl.xlabel("Recall (TP/TP+FN)")
pl.ylabel("Precision (TP/TP+FP)")
pl.grid(True)
pl.show()
    


   