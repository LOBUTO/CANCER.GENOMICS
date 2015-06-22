##Process 1000GENOME Data into table (Mutation frequency in population)
#110914
#Uses data from vcf as in ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf 
#Used language from 150.py
#NOTE!!!:
#    This is to run ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf, however it will run almost completely till it reaches chromosome Y,
#    then it will produce an error, make sure file saves
#    To run chromosome Y make sure to run ALL.chrY.phase3_integrated.20130502.genotypes.vcf separatedly

import vcf 

#Open vcf
FILE_IN1="DATABASES/CANCER_DATA/1000GENOME/2013/ALL.chrY.phase3_integrated.20130502.genotypes.vcf"
VCF_IN=vcf.Reader(open(FILE_IN1)) 

#Write to files. One for structural and the other for snps
FILE_OUT1=open("PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.Y.SNP.CNV", "w")

#First to SNP
for header in ["Chrom","Position","TYPE","REF.ALL","ALT.ALL","EAS.AF","AMR.AF","AFR.AF","EUR.AF","SAS.AF"]:
    FILE_OUT1.write(header+"\t")
FILE_OUT1.write("AF")

for record in VCF_IN:
    print record
    print record.INFO
    
    CHROM=str(record.CHROM)
    POS=str(record.POS)
    REF=str(record.REF)
    
    if "SVTYPE" not in record.INFO:
        TYPE="SNP"
    else:
        TYPE=str(record.INFO["SVTYPE"])
    
    for alt in range(len(record.ALT)):
        EAS=str(record.INFO["EAS_AF"][alt])
        AMR=str(record.INFO["AMR_AF"][alt])
        AFR=str(record.INFO["AFR_AF"][alt])
        EUR=str(record.INFO["EUR_AF"][alt])
        SAS=str(record.INFO["SAS_AF"][alt])
        
        print record.ALT[alt], type(record.ALT[alt]), type(str(record.ALT[alt])), str(record.ALT[alt])
        
        FILE_OUT1.write("\n"+CHROM+"\t"+POS+"\t"+TYPE+"\t"+REF+"\t"+str(record.ALT[alt])+"\t"+
                        EAS+"\t" + AMR +"\t" + AFR +"\t" + EUR +"\t" + SAS +"\t" +
                        str(record.INFO["AF"][alt]))

FILE_OUT1.close()


