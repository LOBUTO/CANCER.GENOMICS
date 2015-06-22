#miRNA project 8
#031814
#Count seed mutation position in mature mirna for cancers

def mirna_mature_position(FILE_IN, FILE_OUT): #of the form 031614_BRCAmaf.mature.mutations
    FILE_IN1=open(FILE_IN)
    MIRNA=[X.split("\t") for X in FILE_IN1.read().splitlines()[1:]]
    FILE_IN1.close()
    
    FILE_OUT1=open(FILE_OUT, "w")
    FILE_OUT1.write("Hugo_Symbol"+"\t"+"Reference_Allele"+"\t"+"Tumor_Seq_Allele2"+"\t"+"Tumor_Sample_Barcode"+"\t"+
                    "Accession"+"\t"+"Name"+"\t"+"Position"+"\t"+"NT")
    
    for record in MIRNA:
        MUT_SEQ=record[14]
        
        #Get mutation position
        X=0
        for nt in MUT_SEQ:
            X=X+1
            if nt.islower():
                POSITION=X
                NT=nt
        
        FILE_OUT1.write("\n"+record[0]+"\t"+record[4]+"\t"+record[5]+"\t"+record[6]+"\t"+record[11]+"\t"+record[13]+
                        "\t"+str(POSITION)+"\t"+NT)
    
    FILE_OUT1.close()

#Run on cancers that have a high population of mature mapped
CANCERS=["BRCA","HNSC","KIRC","LUAD","SKCM","UCEC"]
for cancer in CANCERS:
    mirna_mature_position("DATABASES/CANCER_DATA/miRNA/031614_%smaf.mature.mutations"%cancer,
                          "DATABASES/CANCER_DATA/miRNA/031814_%s_SEED_POS"%cancer)
