#MIXING /CLUSTALO/M1 and /CLUSTALO/PDAU FILES BY ID_X for ALIGNMENT TO M1_X_PDAU

import subprocess

M1=subprocess.check_output("ls",cwd="CLUSTALO/M1/").splitlines() #Get file names from folder
print M1

PDAU=subprocess.check_output("ls",cwd="CLUSTALO/PDAU").splitlines() #Get file names from folder
print PDAU

for i in M1:
    for j in PDAU:
        if str(i)[0:6]==str(j):
            M1A=open("CLUSTALO/M1/%s"%i)
            M1B=M1A.readlines()
            PDAUA=open("CLUSTALO/PDAU/%s"%j)
            PDAUB=PDAUA.readlines()
            SOURCE=M1B[0]+PDAUB[0]
            OUTPUT=open("CLUSTALO/M1_X_PDAU/%s" %i,"w")
            OUTPUT.writelines(SOURCE) #GOT 641, IT MATCHES TO /CLUSTALO/M1 FOLDER 
            OUTPUT.close()
            M1A.close()
            PDAUA.close()
            
            