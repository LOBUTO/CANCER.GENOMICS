#!/bin/bash
#TCGA.CNV.MATRIX.CONSTRUCT.sh
#021615
#Construct a patient x genes matrix using the de-noised data from TCGA and region-based files from processed cnv-matrix
#NOTE: We are using genome hg19 and assuming that those files are in the TCGA cnv.folder
#NOTE: Do in cluster, ANNOVAR mapping may take a very long time to run

#$@ = cnv.folder, annovar directory, cnv.sample.map, outfile
#Remember to initialize in scripts folder

#####Prep CNV files to be mapped######
scripts_folder=$(pwd)
Rscript Function.MAP.CNV.R $1 temp.avinput
mv temp.avinput $2
echo "Done prepping cnv file"
printf "\n"

######Apply annotation througn ANNOVAR######
#Keep in mind that we are running  genome hg19
cd $2 #change to annovar directory to execute perl script
echo "Running ANNOVAR on prepped file - This make take some time"
perl annotate_variation.pl -geneanno -buildver hg19 temp.avinput humandb/
printf "\n"
mv temp.avinput.variant_function $scripts_folder
rm temp.avinput.exonic_variant_function
rm temp.avinput.log 
rm temp.avinput

######Clean up ANNOVAR annotation and Construct Matrix#######
cd $scripts_folder
echo "Cleaning up annovar table and Constructing matrix - This make take some time"
printf "\n"
Rscript Function.ASSIGN.CNV.R temp.avinput.variant_function $3 $1 $4
rm temp.avinput.variant_function
