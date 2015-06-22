#MAP PHASTCON SCORES TO MAF POSITIONS - NOT DOABLE!! TOO LONG TIME - DONE IN BASH, REFER TO COMMENTS IN BOTTOM
#122014
#NOTE!: Script only works if phastcon file is ordered!!!

def command_run (maf_file, phast_file, out_file):
    
    skip_lines=0
    with open(maf_file) as infile:
        next(infile) #Skip header
        
        for line in infile:
            
            #Process each line
            record=line.split("\t")
            chrom=record[0]
            pos=record[1]
            
            #Search for phastcons in phast file
            with open(phast_file) as infile_2:
                
                #Skip previously read lines
                for _ in xrange(skip_lines):
                    next(infile_2)
                
                #Process phast lines
                for phast in infile_2:
                    phast_record=phast.split("\t")
                    phast_chrom=phast_record[0]
                    phast_pos=phast_record[1]
                    
                    #Do while it makes sense to look for it (pos is less than phastcon pos)
                    if (chrom==phast_chrom and int(phast_pos)>int(pos)):
                        break
                    else:
                        skip_lines=skip_lines+1
                    
                    #If match then write with phastcon score
                    if (chrom==phast_chrom and int(pos)==int(phast_pos)):
                        for col in record[:-1]:
                            out_file.write(col+"\t")
                        out_file.write(phast_record[2].rstrip("\n")+"\t")
                        out_file.write(record[-1])
                        break

if __name__ == "__main__":
    
    ####Load Trimer MAF file
    TRIMER_MAF="PIPELINES/METABOLIC.DRIVERS/TABLES/122014.THOUSAND.SNP.TRIMERS.PLUS"
    
    ####Load Phastcons file
    PHAST="DATABASES/PHASTCONS/112714.CHRM.FILTERED.EXONS.100"
    
    ####Load OUTPUT file
    FILE_OUT=open("PIPELINES/METABOLIC.DRIVERS/TABLES/122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST","w")
    FILE_OUT.write("Chrom" +"\t" + "Position" +"\t" + "REF" + "\t" + "ALT" + "\t" + "TRIMER" + "\t" + "MT" +"\t" +"REGION"+"\t"+
                    "PHAST" +"\t"+"AF" +"\n")
    print "pass"
    ####Execute script
    command_run(TRIMER_MAF, PHAST, FILE_OUT)
    
    ####Close file
    FILE_OUT.close()

#awk  'BEGIN {OFS="\t"};{print $1"_"$2, $3, $4, $5, $6, $7, $8}' 122014.THOUSAND.SNP.TRIMERS.PLUS > temp.1
#sort temp.1 > temp.1.sorted

#awk -v OFS="\t" '{print $1 "_" $2, $3}' 112714.CHRM.FILTERED.EXONS.100  > temp.2
#sort temp.2 > temp.2.sorted

#join temp.1.sorted /Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/PHASTCONS/temp.2.sorted > temp.3

#awk -F "_" '{print $1, $2}' temp.3 | awk 'BEGIN {OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > temp.4
#mv temp.4 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST


