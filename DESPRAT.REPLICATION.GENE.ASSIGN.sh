#!/bin/bash
#DESPRAT.REPLICATION.GENE.ASSIGN.sh
#010815
#Assignment of hugo symbols to replication times from the Desprat et al. papter ("Predictable dynamic program of timing of DNA replication in human cells")
#$@ = desprat file, annovar directory, outfile
#despart.file as in 25MS_wig.txt found in the Supplementary_Data_090109 folder
#IMPORTANT: Version of genome is hg18 (36)

######Prep Desprat wigfix file######
scripts_folder=$(pwd)
python desprat.wigfix.prep.py $1 desprat.rt.avinput
mv desprat.rt.avinput $2
echo "Done prepping desprat wigfix file"
printf "\n"

######Apply annotation througn ANNOVAR######
#Again, keep in mind that we are running by default on genome hg18
cd $2 #change to annovar directory to execute perl script
echo "Running ANNOVAR on prepped file"
perl annotate_variation.pl -geneanno -buildver hg18 desprat.rt.avinput humandb/
printf "\n"
mv desprat.rt.avinput.variant_function $scripts_folder
rm desprat.rt.avinput.log 
rm desprat.rt.avinput
rm desprat.rt.avinput.exonic_variant_function

######Clean up annovar annotation##########
echo "Cleaning up annovar table"
cd $scripts_folder
Rscript desprat.annovar.cleanup.R desprat.rt.avinput.variant_function $3
rm desprat.rt.avinput.variant_function

