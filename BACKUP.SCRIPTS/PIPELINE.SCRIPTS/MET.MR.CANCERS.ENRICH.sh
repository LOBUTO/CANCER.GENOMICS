#!/bin/bash
#MET.MR.CANCERS.ENRICH.sh
#012015

AML="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/AML/Somatic_Mutations/WUSM__IlluminaGA_DNASeq/Level_2/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#BRCA="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
#COAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/COAD/100914/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#GBM="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/GBM/100814/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#LUAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/LUAD/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#LUSC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/LUSC/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
HNSC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/HNSC/8998ada8-b0cb-4cf2-bdd0-7321f8904f13/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
KIRC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/KIRC/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq_automated/Level_2/broad.mit.edu__IlluminaGA_automated_DNA_sequencing_level2.maf"
OV="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/OV/101014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#PRAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/PRAD/970b92f0-04d7-4021-980c-9a0c32ecadf3/Somatic_Mutations/BI__IlluminaGA_DNASeq_curated/Level_2/broad.mit.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
READ="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/READ/101014/Somatic_Mutations/BCM__SOLiD_DNASeq/Level_2/hgsc.bcm.edu__ABI_SOLiD_DNA_Sequencing_level2.maf"
SKCM="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/SKCM/608af7c2-400c-4102-bb99-eaced152f04f/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
#STAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/STAD/broad.mit.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
#UCEC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/UCEC/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"

for cancer in $OV
do
    #Obtain cancer type
    TYPE=$(echo $cancer | awk -F "/" '{print $13}')
    echo $TYPE

    #Folder to store
    FILE_OUT="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/$TYPE/012015.$TYPE.MET.HUGO.ENRICHMENT"

    #Run command
    Rscript /Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/SCRIPTS/Function.MAF.METABOLITE.MR.R \
    $cancer \
    /Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/TABLES/011915.EXON.COORDIANTES.FIXED.OVERLAPS \
    /Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds \
    $FILE_OUT
done