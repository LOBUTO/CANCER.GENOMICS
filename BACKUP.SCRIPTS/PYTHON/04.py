#Format i.ligdist5 from capra folder

NAME="1mg5_A"

x="capra_bioinf08_data/1.1.1.1_group0_PF00106_30-85~1.1.1.35_group0_PF00106_30-85"
x=x+".ligdist5" #we could append extension to variable filename

LD5=open(x)
LD5A=[]
for i in LD5:
    LD5A.append(i.split()) #i.split() produces a list of the things being splitted, append in this case is appending lists

print LD5A

LD5B=[]
for i in LD5A:
    for j in i[0:1]:
        if (j!="#" and j!="*"):
            LD5B.append(i) #get only ID_X + all residues within 5A of ligand (with respect to i.aln starting at 0)
        

for i in LD5B: print i

for i in LD5B:
    if i[0]==NAME:
        W=i[1:2] #get first residue number, THIS IS THE ONLY NUMBER I NEED, BUT NEED TO CONNECT TO SPECIFIC "NAME"

print W