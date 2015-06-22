#Reads 1000G vcf files and builds frequency table per mutation at the gene level

import vcf 

#Open 1000G files
FILES_IN=["DATABASES/CANCER_DATA/1000GENOME/ALL.wex.illumina_exome_svm_consensus_v1.20110521.snps.exome.sites.vcf",
          "DATABASES/CANCER_DATA/1000GENOME/ALL.wex.svm_consensus_SOLID_306_hardFiltered.20110521.snps.exome.sites.vcf"]

#Create dict to store
DICT={}

#Work on files
for FILE in FILES_IN:
    
    #Load file
    VCF_IN1=vcf.Reader(open(FILE))
    
    #Keep in mind that one bp mutation can cause a different amino acid change in multiple isoforms of the same gene
    for record in VCF_IN1:
        if "exon" in record.INFO["ANNO"][0]:
            ANNO_LIST=record.INFO["ANNO"]
            ANNO_LIST=filter(lambda x: len(x)>0, ANNO_LIST) #Removes empty strings in list 
            ANNO_LIST_FIRST=ANNO_LIST[0] #Gets only first item since is the only one we need
            
            #Get table elements
            HUGO=str(ANNO_LIST_FIRST.split(":")[1])
            CHROM=str(record.CHROM)
            POS=str(record.POS)
            REF=str(record.REF)
            TOTAL_ALLELES=record.INFO["AN"]
            CLASS=str(ANNO_LIST_FIRST.split(":")[0])
            
            #Write entry in dict for each alternate allele
            ALTERNATE_COUNT=len(record.ALT)
            for alt in range(ALTERNATE_COUNT):
                
                ALT=str(record.ALT[alt])
                ALT_ALLELES=record.INFO["AC"][alt]
                
                #Get combined key
                COMBINED_KEY="$".join([HUGO, CHROM, POS, REF, ALT, CLASS])
                
                #Store or create depending on presence
                if COMBINED_KEY in DICT:
                    DICT[COMBINED_KEY]["ALT_ALLELES"]=DICT[COMBINED_KEY]["ALT_ALLELES"] + ALT_ALLELES 
                    DICT[COMBINED_KEY]["TOTAL_ALLELES"]=DICT[COMBINED_KEY]["TOTAL_ALLELES"] + TOTAL_ALLELES
                
                else:
                    DICT[COMBINED_KEY]={}
                    DICT[COMBINED_KEY]["ALT_ALLELES"]=ALT_ALLELES
                    DICT[COMBINED_KEY]["TOTAL_ALLELES"]=TOTAL_ALLELES
                    
                #print HUGO, CHROM, POS, REF, record.ALT[alt], TOTAL_ALLELES, record.INFO["AC"][alt], CLASS
                
#Open files to write to
FILE_OUT1=open("DATABASES/CANCER_DATA/1000GENOME/062714_1000G_EXON", "w")
FILE_OUT1.write("Hugo_Symbol"+"\t"+"B.Chrom"+"\t"+"B.Position"+"\t"+"REF.BASE"+"\t"+ "ALT.BASE"+"\t"+ 
                "TOTAL.ALLELES"+"\t"+ "ALT.ALLELES"+"\t"+"Variant_Class")     

#Write files
for entry in DICT.keys():
    ENTRY=entry.split("$")
    
    FILE_OUT1.write("\n"+ENTRY[0]+"\t"+ENTRY[1]+"\t"+ENTRY[2]+"\t"+ENTRY[3]+"\t"+ENTRY[4]+"\t"+
                    str(DICT[entry]["TOTAL_ALLELES"])+"\t"+str(DICT[entry]["ALT_ALLELES"])+"\t"+ENTRY[5])
    
FILE_OUT1.close()
