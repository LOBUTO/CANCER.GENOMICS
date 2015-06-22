#USING R FROM PYTHON

from rpy2 import robjects #For regular r functions
from rpy2.robjects.packages import importr #For dev_off()

seq=robjects.r['seq']
mean=robjects.r['mean']
plot=robjects.r['plot']
png=robjects.r['png']
rnorm=robjects.r['rnorm']
barplot=robjects.r['barplot']
length=robjects.r['length']
table=robjects.r['table']
c=robjects.r['c']
grdevices=importr('grDevices')

X=seq(1,5)
print X

Y=seq(0,2,length=5)

grdevices.png("rplot1.png")
A=barplot(Y)
grdevices.dev_off() #to call off graphing mode to be able to graph next plot

A=c('a','b','d','e','b')
print table(A)

barplot(table(A),density=4, ylab="frequency", xlab="amino acid", space=0.5, ylim=c(0,2))

barplot(table(A))


raw_input()