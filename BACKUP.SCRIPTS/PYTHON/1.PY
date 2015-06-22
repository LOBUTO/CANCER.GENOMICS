#GET ALL CAPRA ID_EC with redundancy
import os
print os.getcwd() #this and previous are used to find out what is the working directoty

import glob #function to find all pathnames matching a specified pattern
aln=[]
D=[]
for name in glob.glob("capra_bioinf08_data/*.aln"): #Return a list of pathnames that match pathname specified
    aln.append(name) #creates list aln of filenames *.aln strings

for i in aln:
    data=open(i) #binary to look for flags
    array=[]
    for i in data:  #make into a readable list of elements
        array.append(i)
    A=[]
    for i in array:
        A.append(i.split(" ")) #separate each element into a list so we create a matrix, str.split() takes a string and makes into a list, that's why append works, append(list) takes lists only
    B=[]
    for j in A:
        B.append(j[0:1]) #after checking which columns have it, make those into a new matrix
    C=[]
    for i in B:
        for j in i: #because you cannot split the element but the string inside the element
            C.append(j.split("|")) #separate the string into elements per vector and make a new matrix out of them
   
    for i in C:
        D.append(i[0:1]+i[2:3]) #join the ones you want, protein ID and E.C.

for i in D:
    print i
    
print len(D)