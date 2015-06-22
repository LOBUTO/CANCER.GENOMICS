#1000G.TRIMER.PHASTCON.ASSIGN.sh
#010715
#$@ = trimers.plus.file, phastcon.file, outfile
#trimers file i.e. 122014.THOUSAND.SNP.TRIMERS.PLUS
#phastcon file i.e. 112714.CHRM.FILTERED.EXONS.100

awk  'BEGIN {OFS="\t"};{print $1"_"$2, $3, $4, $5, $6, $7, $8}' $1 > temp.1
sort temp.1 -o temp.1
echo "Done sorting trimers"
printf "\n"

awk -v OFS="\t" '{print $1 "_" $2, $3}' $2  > temp.2
sort temp.2 -o temp.2
echo "Done processing phastcon"
printf "\n"

echo "Joining files"
join temp.1 temp.2 > temp.3
rm temp.1
rm temp.2
echo "Done joining files"
printf "\n"

awk -F "_" '{print $1, $2}' temp.3 | awk 'BEGIN {OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > $3
rm temp.3
echo "Done"

#Notes:
#On 010715 got message: "./1000G.TRIMER.PHASTCON.ASSIGN.sh: line 19: 3: command not found", but did not find error on output