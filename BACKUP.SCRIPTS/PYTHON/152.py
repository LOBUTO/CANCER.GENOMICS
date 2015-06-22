    #Function to process Phastcons

import subprocess

FOLDER="DATABASES/PHASTCONS/100.TRACK"

FILES=subprocess.check_output("ls", cwd=FOLDER).splitlines()
FILES=filter(lambda x: "file" in x, FILES)

FILE_OUT=open("DATABASES/PHASTCONS/070914_PHASTCONS_PROCESSED","w")
FILE_OUT.write("CHROM"+"\t"+"POSITION"+"\t"+"PHASCONS")

for record in FILES:
    print record
    
    #Load file
    FILE_TEMP=open(FOLDER+"/"+record)
    PHAST_FILE=FILE_TEMP.read().splitlines()
    FILE_TEMP.close()
    
    #Process file
    CHROMOSOME=PHAST_FILE[0].split()[1].split("=")[1]
    START=int(PHAST_FILE[0].split()[2].split("=")[1])
    STEP=int(PHAST_FILE[0].split()[3].split("=")[1])
    
    COUNT=0
    for line in PHAST_FILE[1:]:
        FILE_OUT.write("\n"+CHROMOSOME+"\t"+str(START+COUNT)+"\t"+line)
        
        #Update count with step size
        COUNT=COUNT+STEP

FILE_OUT.close()