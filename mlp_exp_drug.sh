#!/bin/bash
#mlp_exp_drug.sh

script="cgp_exp_drugs_mlp.py"
script="cgp_exp_drugsbatch_mlp.py"
tf="True"

for keep_prob in 0.5 0.9 #0.3 0.5 0.9
do
	for e in 2
	do
		for g in 10
		do
			for d in 8_16
			do
				for arch in 2_4 #1_2 1_2_4
				do
					for drug in Erlotinib ZSTK474 CUDC-101 Belinostat Sunitinib CAY10603 \
						TPCA-1 OSI-027 BX-912 BX-795 XL-880 CAL-101 Vorinostat \
						17-AAG 681640 A-443654 A-770041 ABT-263 ABT-869 ABT-888 \
						AC220 AG-014699 AICAR AMG-706 AP-24534 AR-42 AS601245 \
						AS605240 AT-7519 ATRA AUY922 AV-951 AZ628 AZD-0530 AZD6244 \
						AZD6482 AZD7762 AZD8055 Afatinib Axitinib BI-2536 BIX02189 \
						BMN-673 BMS-345541 BMS-509744 BMS-536924 BMS-708163 BMS-754807 \
						Bexarotene Bicalutamide Bosutinib CCT007093 CCT018159 \
						CEP-701 CGP-082996 CGP-60474 CH5424802 CHIR-99021 CI-1040 \
						CMK CP466722 CP724714 Camptothecin Cetuximab Crizotinib \
						Cyclopamine DMOG Dabrafenib Dasatinib EKB-569 EX-527 Elesclomol \
						Embelin Etoposide FK866 FMK FR-180204 FTI-277 GDC0941 GNF-2 \
						GSK-1904529A GSK1070916 GSK2126458 GSK269962A GSK429286A GSK690693 \
						GW-2580 GW843682X Gefitinib HG-6-64-1 I-BET-762
					do
						for slr in 0.00001 # 0.1 0.01 0.001 0.0001 0.00001
						do
							for et in True False #0.1 0.001 0.0001 0.00001 0.000001
							do
								script_name="/tigress/zamalloa/GIT/${script}"
								file_name="${e} ${g} ${d} ${keep_prob} ${arch} ${drug} ${slr} ${tf} ${et}"

								export script_name file_name
								sbatch /tigress/zamalloa/GIT/gs_encode.cmd
							done
						done
					done
				done
			done
		done
	done
done