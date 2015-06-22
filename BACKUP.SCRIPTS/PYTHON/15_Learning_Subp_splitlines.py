#SAVE!! i.splitlines()

import subprocess


A=subprocess.check_output("ls",cwd="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/CLUSTALO/M1/").splitlines()
print len(A)

print 5

B=subprocess.check_output("ls",cwd="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/CLUSTALO/PDAU/").splitlines()
print len(B)

C=[]
for i in A:
    if all(i!=j for j in B):
        C.append(i)

print 6
print C    