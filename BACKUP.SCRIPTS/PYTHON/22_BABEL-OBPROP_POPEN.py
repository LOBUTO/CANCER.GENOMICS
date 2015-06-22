#Testing BABEL - OBPROP

import subprocess

A=open("Components-smiles-oe.smi")
B=A.readlines()
for h in B:
    h=h.split()
    if h[1]=="ILE":
        C=" ".join(h)
print C

CC=open("BABEL/TEST1.smi","w")
CC.writelines(C)
CC.close()

D=subprocess.Popen(["/usr/local/bin/obprop", "BABEL/TEST1.smi"],stdout=subprocess.PIPE) #subprocess.PIPE indicates that a pipe to
                                                                                        #the standard stream should be opened
E=D.communicate()[0].splitlines()  #communicate reads data from stdout until end of file is reached
                                    #[0] tells it to make a list of characters
                                    #so we use splitlines() to produce instead a list of strings(lines) from output

for i in E: print i
print E

for i in E:
    if i[0:4]=="logP":
        F=i.split()[1]
    elif i[0:10]=="mol_weight":
        G=i.split()[1]
print F
print G



