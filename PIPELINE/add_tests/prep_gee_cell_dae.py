# prep_gee_cell_dae.py
# Pre-process Geeleher expression data using DAE after batch normalization
#    by GDSC expression data

import sys
sys.path.insert(0, '/home/ubuntu/GIT/PIPELINE')
from main_pipeline import *
from combat_copy import *
from keras_gsea_dae import *

in_folder   = "/home/ubuntu/"
out_folder  = "/home/ubuntu/TEMP_FEATURES/"

# Load, batch normalize and train data for:
for init_decomp in [1,5]:
	# for arch in ["8", "8_8"]:
	for arch in ["2", "2_2", "2_4_8", "4", "4_4", "4_8"]:
		# DOCETAXEL
		cell_exp_train, cell_exp_test = process_cell_data("GEELEHER/DOCETAXEL/exp_data.txt", in_folder)
		cell_train, cell_test, cell_dae_log = cell_features(cell_exp_train, cell_exp_test, in_folder,
		                                                    CGC=True, test_var=None, 
		                                                    gsea_choice="c2setcovertanimoto", 
		                                                    init_decomp=init_decomp, gsea_dae_arch=arch,
		                                                    gsea_self_valid=False)
		with open("{}/doc_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(out_folder, init_decomp, arch), "wb") as handle:
		    pickle.dump([cell_train, cell_test, cell_dae_log], handle)

		# BORTEZOMIB (A)
		cell_exp_train, cell_exp_test = process_cell_data("GEELEHER/BORTEZOMIB/exp_table_a", in_folder)
		cell_train, cell_test, cell_dae_log = cell_features(cell_exp_train, cell_exp_test, in_folder,
		                                                    CGC=True, test_var=None, 
		                                                    gsea_choice="c2setcovertanimoto", 
		                                                    init_decomp=init_decomp, gsea_dae_arch=arch,
		                                                    gsea_self_valid=False)
		with open("{}/bora_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(out_folder, init_decomp, arch), "wb") as handle:
		    pickle.dump([cell_train, cell_test, cell_dae_log], handle)

		# BORTEZOMIB (B)
		cell_exp_train, cell_exp_test = process_cell_data("GEELEHER/BORTEZOMIB/exp_table_b", in_folder)
		cell_train, cell_test, cell_dae_log = cell_features(cell_exp_train, cell_exp_test, in_folder,
		                                                    CGC=True, test_var=None, 
		                                                    gsea_choice="c2setcovertanimoto", 
		                                                    init_decomp=init_decomp, gsea_dae_arch=arch,
		                                                    gsea_self_valid=False)
		with open("{}/borb_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(out_folder, init_decomp, arch), "wb") as handle:
		    pickle.dump([cell_train, cell_test, cell_dae_log], handle)

		# CISPLATIN
		cell_exp_train, cell_exp_test = process_cell_data("GEELEHER/CISPLATIN/exp_table", in_folder)
		cell_train, cell_test, cell_dae_log = cell_features(cell_exp_train, cell_exp_test, in_folder,
		                                                    CGC=True, test_var=None, 
		                                                    gsea_choice="c2setcovertanimoto", 
		                                                    init_decomp=init_decomp, gsea_dae_arch=arch,
		                                                    gsea_self_valid=False)
		with open("{}/cisp_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(out_folder, init_decomp, arch), "wb") as handle:
		    pickle.dump([cell_train, cell_test, cell_dae_log], handle)

print("Done prepping all gee cell features")