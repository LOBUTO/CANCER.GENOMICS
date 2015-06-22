#LEARNING TO PLOT 

import numpy as np #FOR IMPORTING ARRAYS
import pylab as pl #FOR PLOTTING
import scipy.stats #FOR CALCULATIONS SPEARMAN CORRELATIONS

X=[1,2,3,4,5]
Y=[]
for i in X:
    Y.append(i+2)
print Y
"""
pl.plot(X,Y,"r--")#FOR PLOT SPECIFY VARIABLES FIRST (X,Y), THEN YOU COULD ADD LINE STYLE (--) OR COLOR (r)
                    #MAKE SURE TO PUT THEM AS A SINGLE STRING FOLLOWING THE VARIABLES
pl.show() #SHOW PLOT IN SCREEN, CAN SHOW SEQUENTIAL PLOTS

pl.plot(X,Y,"bo")#IN THIS CASE WE CHOOSE BLUE(b) AND SCATTER PLOT CIRCLES (o)
pl.show()

pl.plot(X,Y,"g*")#IN THIS CASE WE CHOOSE GREEN(g) AND SCATTER PLOT ASTERISKS(*)
pl.show()

pl.plot(X,Y,"bs")#IN THIS CASE WE CHOOSE BLUE(b) AND SQUARES (s)
pl.xlabel("Xs") #LABELS FOR AXIS
pl.ylabel("Ys")
pl.title("TESTING GRAPHING") #TITLE OF GRAPH
pl.xlim(X[0]-1,X[-1]+1) #SET RANGE
pl.ylim(Y[0]-1,Y[-1]+1)
pl.show()

pl.plot(X,Y,"g*")#IN THIS CASE WE CHOOSE GREEN(g) AND SCATTER PLOT ASTERISKS(*)
pl.show()
"""

#LOADING NUMERICAL DATA WITH NUMPY AND CALCULATING PEARSON AND SPEARMAN
data=np.loadtxt("SMILES/32_LOG_P") #Load file into an array [rows,columns]
print data
print data[:,0] #PRINTS ALL ROWS OF FIRST COLUMN
pl.plot(data[:,0],data[:,1],"ro")
pl.show()

print scipy.stats.pearsonr(data[:,0],data[:,1])
print scipy.stats.spearmanr(data[:,0],data[:,1])

DATA=open("SMILES/32_LOG_P") #HAVE TO OPEN IT REGULARLY TO GET RID OF ZEROES
DATA1=[]
for i in DATA:
    DATA1.append(i.split())
    
DATA1=filter(lambda x:x[0]!="0.0" and x[1]!="0.0",DATA1) #FILTERING OUT ZEROES
X=[]
Y=[]
for i in DATA1:
    X.append(float(i[0]))
    Y.append(float(i[1]))

pl.plot(X,Y,"b*")
pl.xlim(sorted(X)[0]-0.5,sorted(X)[-1]+0.5)
pl.ylim(sorted(Y)[0]-0.5,sorted(Y)[-1]+0.5)
pl.show()

print scipy.stats.pearsonr(X,Y)
print scipy.stats.spearmanr(X,Y)
