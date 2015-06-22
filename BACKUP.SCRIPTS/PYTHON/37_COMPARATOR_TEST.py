#COMPARATOR TEST

import pylab as pl

FILE_IN=open("SMILES/32_LOG_P2_MEAN")
FILE=[]
for i in FILE_IN:FILE.append(i.split())

FILE_NON_ZERO=filter(lambda x: x[0]!="0.0" and x[1]!="0.0", FILE)
print FILE_NON_ZERO
print len(FILE_NON_ZERO)

POS=0
NEU=0
for i in FILE_NON_ZERO:
    if float(i[0])>0 and float(i[1])>0:
        POS=POS+1
    elif float(i[0])>0 and float(i[1])<0:
        NEU=NEU+1
    elif float(i[0])<0 and float(i[1])>0:
        NEU=NEU+1
    elif float(i[0])<0 and float(i[1])<0:
        POS=POS+1
        
print POS,NEU

Y=[NEU,POS]
N=len(Y)
ind=range(N)
FIG=pl.figure()
ax=FIG.add_subplot(1,1,1)
ax.set_xticks(ind)
ax.set_ylabel("Counts")
ax.bar(ind,Y, width=0.5)
group_labels=["Different change", "Same change"]
ax.set_xticklabels(group_labels)
pl.xlim(-1,2.5)
pl.ylim(0,20)
FIG.autofmt_xdate()
pl.show()

print FILE
print len(FILE)

POS=0
NEU=0
for i in FILE:
    if float(i[0])>0 and float(i[1])>0:
        POS=POS+1
    elif float(i[0])>0 and float(i[1])<0:
        NEU=NEU+1
    elif float(i[0])<0 and float(i[1])>0:
        NEU=NEU+1
    elif float(i[0])<0 and float(i[1])<0:
        POS=POS+1
print POS,NEU