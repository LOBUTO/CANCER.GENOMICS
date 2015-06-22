#Function_process1000G.py
#070714
#Reads 1000G vcf files and builds frequency table per mutation at the amino acid level

import vcf, sys

OUTPUT_FILE=sys.argv[1]
VCF_FILES=sys.argv[2:] #1000G files

print OUTPUT_FILE
print VCF_FILES

#Complementary Dict
DICT_COMP={"A":"T","C":"G","T":"A","G":"C"}

#Create dict to store
DICT={}

#Work on files
for FILE in VCF_FILES:
    
    #Load file
    VCF_IN1=vcf.Reader(open(FILE))
    
    #Keep in mind that one bp mutation can cause a different amino acid change in multiple isoforms of the same gene
    for record in VCF_IN1:
        
        if "exon" in record.INFO["ANNO"][0]:
            ANNO_LIST=record.INFO["ANNO"]
            ANNO_LIST=filter(lambda x: len(x)>0, ANNO_LIST) #Removes empty strings in list 
            
            HUGO=str(ANNO_LIST[0].split(":")[1])
            CHROM=str(record.CHROM)
            TOTAL_ALLELES=record.INFO["AN"]
            CLASS=str(ANNO_LIST[0].split(":")[0])
            
            if len(record.ALT)==1:
                
                ALT_ALLELES=record.INFO["AC"][0]
                
                #Extract all possible mutation site pairs (i.e. "T181A:p.S61T")
                MUTATIONS=list(set([mut.split(":c.")[1] for mut in ANNO_LIST]))
                
                #Append each mutation to dictionary:
                for mutation in MUTATIONS:
                    COMBINED_KEY="$".join([HUGO, CHROM, CLASS, mutation])
                
                #Store or create depending on presence
                if COMBINED_KEY in DICT:
                    DICT[COMBINED_KEY]["ALT_ALLELES"]=DICT[COMBINED_KEY]["ALT_ALLELES"] + ALT_ALLELES 
                    DICT[COMBINED_KEY]["TOTAL_ALLELES"]=DICT[COMBINED_KEY]["TOTAL_ALLELES"] + TOTAL_ALLELES
                
                else:
                    DICT[COMBINED_KEY]={}
                    DICT[COMBINED_KEY]["ALT_ALLELES"]=ALT_ALLELES
                    DICT[COMBINED_KEY]["TOTAL_ALLELES"]=TOTAL_ALLELES
                
            else: #For more than one alternate allele
                
                #Need REF to see if it is in complementary or not to know which allele count to call
                REF=str(record.REF)
                
                #Extract base mutations
                PRE_MUTATIONS=list(set([mut.split(":c.")[1] for mut in ANNO_LIST]))
                
                #Check if complementary or not
                POST_MUTATIONS=[]
                for pair in PRE_MUTATIONS:
                    if pair[0]==REF:
                        POST_MUTATIONS=POST_MUTATIONS+[pair]
                    else:#If complementary need to change it so that we can match it to the count in alt alleles
                        comp_pair=pair.replace(pair[0], DICT_COMP[pair[0]],1)
                        comp_pair=comp_pair.split(":p.")[0][:-1]+DICT_COMP[comp_pair.split(":p.")[0][-1]] + ":p." +comp_pair.split(":p.")[1]
                        POST_MUTATIONS=POST_MUTATIONS+[comp_pair]
                        
                #Extract allele count by presence of base in amino acid mutation calling 
                ALTERNATE_COUNT=len(record.ALT)
                for allele in range(ALTERNATE_COUNT):
                    
                    ALT=str(record.ALT[allele])
                    ALT_ALLELES=record.INFO["AC"][allele]
                    
                    #Check that ALT is present in mutations called to count allele frequency, otherwise don't count
                    MUTATIONS=filter(lambda mut: mut.split(":p.")[0][-1]==ALT ,PRE_MUTATIONS)
                    
                    #Append each mutation to dictionary:
                    for mutation in MUTATIONS:
                        COMBINED_KEY="$".join([HUGO, CHROM, CLASS, mutation])
                    
                    #Store or create depending on presence
                    if COMBINED_KEY in DICT:
                        DICT[COMBINED_KEY]["ALT_ALLELES"]=DICT[COMBINED_KEY]["ALT_ALLELES"] + ALT_ALLELES 
                        DICT[COMBINED_KEY]["TOTAL_ALLELES"]=DICT[COMBINED_KEY]["TOTAL_ALLELES"] + TOTAL_ALLELES
                    
                    else:
                        DICT[COMBINED_KEY]={}
                        DICT[COMBINED_KEY]["ALT_ALLELES"]=ALT_ALLELES
                        DICT[COMBINED_KEY]["TOTAL_ALLELES"]=TOTAL_ALLELES

                
#Open files to write to
FILE_OUT1=open(OUTPUT_FILE, "w")
FILE_OUT1.write("Hugo_Symbol"+"\t"+"B.Chrom"+"\t"+"B.Position"+"\t"+"REF.AA"+"\t"+ "ALT.AA"+"\t"+ 
                "TOTAL.ALLELES"+"\t"+ "ALT.ALLELES"+"\t"+"Variant_Class")     

#Filter out 0 alternate allele count
DICT_KEYS=filter(lambda keys: DICT[keys]["ALT_ALLELES"]!=0, DICT.keys())

#Write files
for entry in DICT_KEYS:
    ENTRY=entry.split("$")
    
    FILE_OUT1.write("\n"+ENTRY[0]+"\t"+ENTRY[1]+"\t"+ENTRY[3].split(":p.")[1][1:-1]+"\t"+
                    ENTRY[3].split(":p.")[1][0]+"\t"+ENTRY[3].split(":p.")[1][-1]+"\t"+
                    str(DICT[entry]["TOTAL_ALLELES"])+"\t"+str(DICT[entry]["ALT_ALLELES"])+"\t"+ENTRY[2])
    
FILE_OUT1.close()
