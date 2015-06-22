#Get M1_O_TOTAL
import subprocess

M1A=subprocess.check_output("ls", cwd="capra_bioinf08_data/").splitlines()
print len(M1A)

#Get unique ID_X from all i.csadist5 files
CSD=[]
ALN=[]
LD=[]
for i in M1A:
    if i[-8:]=="csadist5":
        CSD.append(i[:-8])
    elif i[-3:]=="aln":
        ALN.append(i[:-3])
    elif i[-8:]=="ligdist5":
        LD.append(i[:-8])    

print ALN    
print CSD
print LD

M1=[]
for j in CSD:
    for k in ALN:
        if j==k:
            CSDA=open("capra_bioinf08_data/%scsadist5"%j) #append extension to open file
            CSDB=CSDA.readlines()
            CSDC=[]
            for l in CSDB:
                if l[0]!="#" and l[0]!="*" and l[0]!="\n":
                    CSDC.append(l[0:6]) #List of ID_X in i.csadist5 files
                    
            ALNA=open("capra_bioinf08_data/%saln"%j)
            ALNB=ALNA.readlines()
            ALNC=[]
            for m in ALNB: #Get ID_Xs from i.csadist5 in i.aln
                for h in CSDC:
                    if m[0:6]==h:
                        m=m.split()
                        ALNC.append(m[0][0:6]+" "+m[1]) #ID_X plus skipped sequences
            
            for q in CSDC: #For each ID_X in i.csadist5
                ALNE="" #To store sequence
                for r in ALNC:
                    if q==r[0:6]: #if it is equal to ID_X in ALNC
                        ALNE=ALNE+r[7:]
            
                M1.append(q+" " +ALNE+" "+"%ssdpO"%j+"\n")
 
ALNA.close()
    
for o in LD:
    for k in ALN:
        if o==k:
            LDA=open("capra_bioinf08_data/%sligdist5"%o) #append extension to open file
            LDB=LDA.readlines()
            LDC=[]
            for p in LDB:
                if p[0]!="#" and p[0]!="*" and p[0]!="\n":
                    LDC.append(p[0:6]) #ID_X list from i.ligdist5          
            
            ALNA=open("capra_bioinf08_data/%saln"%o)
            ALNB=ALNA.readlines()
            ALNC=[]
            for m in ALNB: #Get ID_Xs from i.csadist5 in i.aln
                for h in LDC:
                    if m[0:6]==h:
                        m=m.split()
                        ALNC.append(m[0][0:6]+" "+m[1]) #ID_X plus skipped sequences         
            print ALNC
            
            for q in LDC: #For each ID_X in i.csadist5
                ALNE="" #To store sequence
                for r in ALNC:
                    if q==r[0:6]: #if it is equal to ID_X in ALNC
                        ALNE=ALNE+r[7:]            

                M1.append(q+" " +ALNE+" "+"%ssdpO"%o+"\n")

ALNA.close()

M1.sort(key=lambda X:X[0:6])
#for i in M1: print i

M1_CSADIST=open("M1_O_TOTAL", "w")
M1_CSADIST.writelines(M1)

M1_CSADIST.close()