#TCGA.TRIMER.ASSIGN.sh
#123014
#$@ = TCGA.maf, file.out, fasta.file, phastcon.file
#fasta.file as in human_g1k_v37.fasta
#phastcon.file as in 112714.CHRM.FILTERED.EXONS.100

########Pre-process TCGA.maf###############
Rscript TCGA.PRE.PROCESS.R $1
printf "\n"

#######Assign trimers to each position#####
python TCGA.TRIMER.SEARCH.py temp.1.csv temp.2.csv $3
rm temp.1.csv
echo "Done processing trimers"
printf "\n"

#####Assign PHASTCON scores################
echo "Sorting files"
printf "\n"

awk  'BEGIN {OFS="\t"};{print $4 "$" $1, $2, $3, $5, $6, $7, $8, $9}' temp.2.csv > temp.3
sort temp.3 -o temp.3
rm temp.2.csv

awk -v OFS="\t" '{print $1 "$" $2, $3}' $4  > temp.4
sort temp.4 -o temp.4
echo "Done sorting files"
printf "\n"

echo "Assigning phastcons scores"
printf "\n"

join temp.3 temp.4 > temp.5
rm temp.3
rm temp.4

awk -F "$" '{print $1, $2}' temp.5 | awk  '{for (i=1;i<=NF;i++){printf "%s\t",  $i};printf "\n"; }' > $2
rm temp.5
echo "Done"

#In 123114 produced 123114.BRCA.MAF.TRIMERS.csv with breast cancer data from TCGA
