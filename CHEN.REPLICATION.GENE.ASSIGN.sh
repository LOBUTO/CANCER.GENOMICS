#!/bin/bash
#CHEN.REPLICATION.GENE.ASSIGN.sh
#010815
#Assignment of hugo symbols to replication times from the Chen et al. paper ("Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes")
#$@ = chen.file, annovar directory, outfile
#chen.file as in S50.100kbp.Non-overlappingWindows.hg18.wig
#IMPORTANT: Version of genome is hg18 (36)
#ALSO: Keep in mind that the chen file does not have replication time for "X" and "Y" chromosomes!!!

######Prep Chen wigfix file######
scripts_folder=$(pwd)
python chen.wigfix.prep.py $1 chen.rt.avinput
mv chen.rt.avinput $2
echo "Done prepping chen wigfix file"
printf "\n"


######Apply annotation througn ANNOVAR######
#Again, keep in mind that we are running by default on genome hg18
cd $2 #change to annovar directory to execute perl script
echo "Running ANNOVAR on prepped file"
perl annotate_variation.pl -geneanno -buildver hg18 chen.rt.avinput humandb/
printf "\n"
mv chen.rt.avinput.variant_function $scripts_folder
rm chen.rt.avinput.log 
rm chen.rt.avinput
rm chen.rt.avinput.exonic_variant_function

######Clean up annovar annotation##########
echo "Cleaning up annovar table"
cd $scripts_folder
Rscript chen.annovar.cleanup.R chen.rt.avinput.variant_function $3
rm chen.rt.avinput.variant_function

