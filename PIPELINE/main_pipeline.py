# main_pipeline.py

import sys, gc
from absl import flags

sys.path.insert(0, '/home/ubuntu/GIT')
from paper_1 import *

sys.path.insert(0, '/home/ubuntu/GIT/PIPELINE')
from combat_copy import *
from keras_gsea_dae import *
from dae_cross import *
from ml_pgm_models import *

def load_train_exp(in_folder):
	# Loads expression data obtained from GDSC as samplesxgenes
	x = pd.read_csv("{}CGP_FILES/070818_cgp_exp.txt".format(in_folder), sep="\t")
	return x

def load_test_exp(cell_array):
	# Load external expression dataset file: samplesxgenes
	x = pd.read_csv(cell_array, sep="\t")
	return x

def exp_consolidate(a_exp, b_exp):
	# Consolidates two expression matrices into singles one aligning common genes
	# Matrices have to be in the form samples x genes

	common_genes = list(set(a_exp.columns.values).intersection(b_exp.columns.values))
	print("Number of genes across datasets: {}".format(len(common_genes)))

	# Fix sample names for duplicity
	a_exp.index  = a_exp.index + "_first"
	b_exp.index  = b_exp.index + "_second"

	main_exp     = pd.concat([a_exp[common_genes], b_exp[common_genes]])

	return main_exp

def batch_list(a_exp, b_exp):
	# Creates batch list pandas frame as sample labels in index and batch correspondence in first column
	# NOTE: Constructed in the order the tables were given
	# NOTE: Assume that names have already been fixed for duplicity (no duplicate sample names across tables)

	batch = pd.DataFrame({"batch": ["first"]*a_exp.shape[0] + ["second"]*b_exp.shape[0]},
						 index=list(a_exp.index) + list(b_exp.index))
	return batch

def post_process_batch_norm(batch_data):
	# Clean up and separate batch normalized processed data

	first  = batch_data.filter(regex="_first$", axis=1)
	second = batch_data.filter(regex="_second$", axis=1)

	first.columns  = [i.rstrip("_first") for i in list(first.columns)]
	second.columns = [i.rstrip("_second") for i in list(second.columns)]

	return first, second

def process_cell_data(cell_array, in_folder):
	# Function to normalize gene expression data from test and train sources
	# Returns normalized expression matrices as samples x genes

	# Load train data expression array
	train_exp  = load_train_exp(in_folder)  # By default we load data from GDSC
	test_exp   = load_test_exp(cell_array)

	main_exp   = exp_consolidate(train_exp, test_exp)

	# Batch normalize expression arrays
	batch_var  = batch_list(train_exp, test_exp)
	batch_norm = combat(main_exp.T, batch_var["batch"])

	# Clean up
	train_exp, test_exp = post_process_batch_norm(batch_norm)

	return train_exp.T, test_exp.T

def change_smiles_morgan_fingerprint(smiles_list):
	# smiles_list: list of smiles strings
	# Uses rdkit
	from rdkit import Chem
	from rdkit.Chem import AllChem

	main_list = []
	for i in smiles_list:
		m1  = Chem.MolFromSmiles(i)
		fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=1024)
		fp1 = fp1.ToBitString()
		main_list.append([int(i) for i in str(fp1)])

	return main_list

def change_smiles_pubchem_fingerprint(smiles_list):
	# smiles_list: list of smiles strings
	# Uses pubchempy
	import pubchempy as pcp

	main_list = []
	for i in smiles_list:
		p_id        = pcp.get_cids(identifier=i, namespace="smiles")
		fingerprint = pcp.Compound.from_cid(p_id).cactvs_fingerprint
		main_list.append([int(i) for i in str(fingerprint)]) 

	return main_list

def change_pubchemid_smiles(id_list):
	# id_list: list of integer ids according to pubchem
	import pubchempy as pcp

	main_list = []
	for i in id_list:
		main_list.append(pcp.Compound.from_cid(i))

	return main_list

def construct_fingerprint_features(test_drug):
	# Constructs feature matrix out of 2-column pandas table of names and smiles

	all_smiles = list(test_drug.smiles)
	all_drugs  = list(test_drug.Compound)
	features   = change_smiles_morgan_fingerprint(all_smiles)

	feat_table = pd.DataFrame(features, index=all_drugs,
							  columns = ["s_%s"%j for j in xrange(len(features[0]))] )

	return feat_table

def construct_mordred_features(table_in):
	# Constructs feature matrix from mordred physico-chemical features
	#  out of 2-column pandas table of names and smiles [Compound, smiles]
	from rdkit import Chem
	from mordred import Calculator, descriptors

	# Create descriptors
	calc       = Calculator(descriptors, ignore_3D=False)

	# Get features
	all_smiles = list(table_in.smiles)
	all_drugs  = list(table_in.Compound)	
	mols       = [Chem.MolFromSmiles(smi) for smi in all_smiles]

	# Clean up
	feat_table = calc.pandas(mols)
	feat_table = feat_table.select_dtypes(["number"])
	feat_table.index = all_drugs

	return feat_table

def load_gdsc_smiles(in_folder):
	# Loads latest gdsc updated smiles data as two-column table
	# Extracts canonical smiles from pubchem_ID - using PUG-REST API!!
	#import pubchempy as pcp
	import urllib2

	# Load files
	file_in    = "{}CGP_FILES/pubchem_id.csv".format(in_folder)
	file_in_updated = "{}CGP_FILES/gdsc_drug_list.csv".format(in_folder)

	x          = pd.read_csv(file_in, header=None, names=["Compound", "count", "ID","N", "DROP"]).drop_duplicates() # SOURCE 1
	x_update   = pd.read_csv(file_in_updated, usecols=[1,5], header=0, names=["Compound", "ID"]).drop_duplicates() # SOURCE 2

	# Clean up
	x["Compound"] = [i.rstrip("'") for i in list(x.Compound)]
	x.ID          = x.ID.astype("str")
	x_update      = x_update.query("ID!='none'").query("ID!='several'")
	x_update      = x_update.dropna()

	# Combine with updated
	x_update   = x_update.loc[~x_update.Compound.isin(list(x.Compound))]
	x          = pd.concat([x[["Compound", "ID"]], x_update])

	#NOTE: Manual clean up of IDs due to spotted technical duplicates#
	remove_ids = ["681640", "53298813", "57370134", "16760671", "447912", "16683866", "6851740"]
	x          = x.loc[~x.ID.isin(remove_ids)]
	x          = pd.concat([x, pd.DataFrame([["T0901317", "447912"]], columns=x.columns)])
	x          = x.drop_duplicates()
	##################################################################

	# Extract smiles using pubchem IDs
	all_ids    = list(x.ID)
	id_string  = ",".join(all_ids)

	url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/CanonicalSMILES/CSV".format(id_string)
	response   = urllib2.urlopen(url)
	response   = pd.read_csv(response, names=["ID", "smiles"])

	main_table = pd.merge(x, response, on="ID")[["Compound", "smiles"]].drop_duplicates()

	# Clean up GDSC names
	main_table["Compound"] = [Function_drug_name_zhang_to_gdsc(i, "gtz") for i in list(main_table.Compound)]
	main_table = main_table.drop_duplicates()

	# Return
	return main_table

def process_drug_smiles(drug_smiles, in_folder):
	# Process smiles-based features for train and test drug compounds
	# - Takes in test drugs smiles and converts them into fingerprints and mordred features
	# - Load internal train drug dataset for fingerprints
	# NOTE: Fingerprints are based on pubchempy, a set of 881 features

	# Load smiles data
	train_smiles   = load_gdsc_smiles(in_folder)
	test_smiles    = pd.read_csv(drug_smiles, sep="\t", header=None, names=["Compound", "smiles"])

	# Obtain fp train data
	# train_fp       = pd.read_csv("{}CGP_FILES/pubchem_smiles_updated.txt".format(in_folder), index_col=0)
	# train_fp.index = [Function_drug_name_zhang_to_gdsc(i, "gtz") for i in list(train_fp.index)]
	# train_fp       = train_fp.rename_axis("Compound")
	train_fp       = construct_fingerprint_features(train_smiles)

	# Extract mordred train data - NOTE: AUTOMATIZED FROM STORED FILE
	with open("{}MORDRED_FILES/gdsc_mordred.pickle".format(in_folder), "rb") as handle:
		train_mordred = pickle.load(handle)
	# train_mordred  = construct_mordred_features(train_smiles)

	# Extract fp and mordred test features
	test_fp        = construct_fingerprint_features(test_smiles)
	test_mordred   = construct_mordred_features(test_smiles)

	# Consolidate mordred features
	common_mordred = list(set(train_mordred.columns).intersection(test_mordred.columns))
	train_mordred  = train_mordred[common_mordred]
	test_mordred   = test_mordred[common_mordred]

	# Return
	return train_fp, train_mordred, test_fp, test_mordred

def process_giant_features(cell_exp_train, cell_exp_test, in_folder, th=0.98):
	# Produces compatible CGC-GIANT expression features

	genes = process_cgc_file("{}CGC_FILES/cancer_gene_census.csv".format(in_folder))[1]

	# Filter GIANT genes by threshold
	giant     = process_giant_entrez("{}GIANT_FILES/all_tissues_top".format(in_folder), 
									 "{}HGNC_FILES/genes_entrez.txt".format(in_folder))
	giant     = giant.loc[giant.prob>=th]

	# Filter by cgc genes
	giant     = giant.loc[(giant.gene_1.isin(genes)) | (giant.gene_2.isin(genes))]
	genes     = list(set( list(giant.gene_1) + list(giant.gene_2)))

	# Apply to expression values and extract relevant genes
	exp_genes = set(list(cell_exp_train.columns)).intersection(list(cell_exp_test.columns))
	exp_genes = list(exp_genes.intersection(genes))

	# Assess compabitilibity
	cell_exp_train = cell_exp_train[exp_genes]
	cell_exp_test  = cell_exp_test[exp_genes]

	# Return
	return cell_exp_train, cell_exp_test

def process_var_features(test_var, in_folder):
	# Produces compatible variant feature matrices based on amino acid variations

	# Process train data
	train_var   = process_gdsc_variant("{}CGP_FILES/gdsc_WES_variants.csv".format(in_folder), 
									output="matrix", th=5)

	# Process test data
	test_var    = pd.read_csv(test_var, sep="\t", index_col=0)

	# Assess compatibility
	common_vars = list(set(train_var.columns).intersection(test_var.columns))

	# Return if common variants found
	if len(common_vars)>2:
		print("{} common variants found between training and test data".format(len(common_vars)))
		train_var, test_var = train_var[common_vars], test_var[common_vars]

		# Add identifier
		train_var.columns   = [i+"_var" for i in list(train_var.columns)]
		test_var.columns    = [i+"_var" for i in list(test_var.columns)]

		return train_var, test_var
	else:
		print("No common variants found")
		return None, None

def cell_features(cell_exp_train, cell_exp_test, in_folder,
				  gsea_choice="c2setcovertanimoto", CGC=False, train_var=None,test_var=None, 
				  init_decomp=5, gsea_dae_arch="2", gsea_self_valid=True):
	# Function to extract ALL types of cell-related features (expression, variation)
	# *_exp_train: as panda matrices as samples x genes
	# test_var: matrix text file of samples x features(gene_proteinmut)
	# 	- protein mutation as "p.A#B", where:
	#		- A: Original amino acid, B: Final amino acid, # amino acid position
	#		- Can place "*" in place of "B" 

	# Obtain GSEA-DAE features
	# NOTE: Track "dae_log" for DAE stats
	train_dae, test_dae, dae_log = gsea_autoencoder(cell_exp_train, cell_exp_test, in_folder,
												gsea=gsea_choice, filter_th=int(10), arch=gsea_dae_arch, labeller="_gs", 
												init_decomp=init_decomp, self_valid=gsea_self_valid)

	# Obtain CGC-GIANT features (if necessary)
	if CGC==True:
		th    = 0.98
		train_cgc, test_cgc = process_giant_features(cell_exp_train, cell_exp_test, in_folder, th)
	else:
		train_cgc, test_cgc = None, None

	# Obtain Variant features (if available)
	if test_var is not None:
		print("Using sample's mutation information as features")
		train_var, test_var = process_var_features(test_var, in_folder)

		# Clean up for name and index matching
		train_var = cell_gdsc_zhang_rename_index(train_var, "gtz")
		train_var = train_var.loc[set(cell_exp_train.index).intersection(train_var.index)]
		test_var  = test_var.loc[set(cell_exp_test.index).intersection(test_var.index)]

	# Consolidate all features
	train_features = pd.concat([train_dae, train_cgc, train_var], axis=1)
	test_features  = pd.concat([test_dae, test_cgc, test_var], axis=1)
	print("Number of cell samples available for training: {}".format(train_features.shape[0]))

	# Clean up and return
	train_features = train_features.fillna(0) # To account for those that do no contain genomic variance at positions (train_var)
	test_features  = test_features.fillna(0) # To account for those that do no contain genomic variance at positions (test_var)

	return train_features, test_features, dae_log

def encoded_target_giant(test_target, in_folder, th=0.7):
	# Function to produce compatible Target-GIANT DAE feature matrices

	# Load train GIANT features
	# NOTE: If this process takes too long we should have the data ready to go
	gdsc_target    = process_drug_target("{}CGP_FILES/Screened_Compounds.csv".format(in_folder))
	giant          = process_giant_entrez("{}GIANT_FILES/all_tissues_top".format(in_folder), 
										  "{}HGNC_FILES/genes_entrez.txt".format(in_folder)) #["gene_1", "gene_2", "prob"]

	train_features = giant_gene_features(giant, gdsc_target, th=th)

	# Load test GIANT features
	test_target    = pd.read_csv(test_target, sep="\t", names=["Compound", "target"])
	test_features  = giant_gene_features(giant, test_target, th=th)

	# Consolidate features
	print("Train drug targets: {}".format(train_features.shape[1]))
	print("Test drug targets: {}".format(test_features.shape[1]))
	features       = list(set(train_features.columns).intersection(list(test_features.columns)))
	print("Common drug target features: {}".format(len(features)))
	train_features = train_features[features]
	test_features  = test_features[features]

	# Obtain Autoencoded features
	# NOTE: Track dae_log GIANT stats
	train_features, test_features, dae_log = cross_feature_dae(train_features, test_features, arch="4_16", labeller="_giant")

	# Return
	return train_features, test_features, dae_log

def drug_features(drug_smiles, in_folder, train_target=None, test_target=None):
	# Function to extract ALL types of drug-related features
	# test_target: two-column tab-separated file of drug-target

	# Extract smiles-related features first
	train_fp, train_mordred, test_fp, test_mordred = process_drug_smiles(drug_smiles, in_folder)

	# Autoencode mordred-based features
	train_mordred, test_mordred, mordred_dae_log   = cross_feature_dae(train_mordred, test_mordred, arch="8_16", labeller="_mordred")

	# Obtain Target-GIANT DAE features (if available)
	giant_dae_log     = None
	if test_target is not None:
		print("Using Drug-target information")
		th = 0.7
		train_target, test_target, giant_dae_log   = encoded_target_giant(test_target, in_folder, th=th)
		giant_dae_log = giant_dae_log.assign(source="Drug_Giant_dae")

		# Fix for common names and index matching
		train_target = drug_gdsc_zhang_rename_index(train_target, "gtz")
		train_target = train_target.loc[set(train_fp.index).intersection(train_target.index)]
		test_target  = test_target.loc[set(test_fp.index).intersection(test_target.index)]

	# Consolidate all features
	train_features = pd.concat([train_fp, train_mordred, train_target], axis=1)
	test_features  = pd.concat([test_fp, test_mordred, test_target], axis=1)

	# Clean up and return
	train_features = train_features.fillna(0) # To account for those that do no contain drug target information (train_target)
	print("Number of drug samples available for training: {}".format(train_features.shape[0]))
	test_features  = test_features.fillna(0) # To account for those that do no contain drug target information (test_target)
	log_table      = pd.concat([mordred_dae_log.assign(source="Drug_mordred_dae"),
								giant_dae_log])

	return train_features, test_features, log_table

def pipeline(in_folder, cell_array, drug_smiles, target, drug_target=None, sample_mut=None):
	# cell_array: text array of samples x genes
	# drug_smiles: 2 column tab-separated file of compounds and smiles, tab-separated
	#	ie. CTRP_FILES/drug_smiles.txt
	# target: 2 column tab-separated file of compounds-cell_name combinations to predict

	# Batch normalize gene expression data - samples x genes
	cell_exp_train, cell_exp_test       = process_cell_data(cell_array, in_folder)

	# Extract cell features
	cell_train, cell_test, cell_dae_log = cell_features(cell_exp_train, cell_exp_test, in_folder,
														CGC=True, test_var=sample_mut, 
														init_decomp=5, gsea_dae_arch="2",
														gsea_self_valid=True)

	# Extract drug features
	drug_train, drug_test, drug_dae_log = drug_features(drug_smiles, in_folder, test_target=drug_target)
	
	# Extract training target variable
	target_table                        = cgp_act_post_process(pd.read_csv("{}CGP_FILES/v17.3_fitted_dose_response.csv".format(in_folder)),
											zscoring=True)

	# Clean up for correct index names
	print("Number of genomic features: {}".format(cell_train.shape[1]))
	print("Number of drug molecular features: {}".format(drug_train.shape[1]))
	target_table, cell_train, drug_train = drug_exp_lobico_parse_v3(target_table, cell_train, drug_train, "gtz")

	# Model train data and predict
	test_prediction, deep_pgms_log       = deep_pgms(target_table, cell_train, drug_train,
													cell_test, drug_test, target,
													arch="32_32_32_32", l2_reg=0.0001, keepprob=0.5)

	# Return prediction and logs
	return test_prediction, [cell_dae_log, drug_dae_log, deep_pgms_log]

def write_output(out_file, output, log):
	# Write main predictions
	output.to_csv(out_file, sep="\t", index=False)

	# Write log files
	with open("{}.log.pickle".format(out_file), "wb") as handle:
		pickle.dump(log, handle)

def run():

	# Arguments
    cell_array  = sys.argv[1] # cell x genes file, such as ctrp_exp_2.1.txt
    drug_smiles = sys.argv[2] # 2 column file of drug_name - smiles
    target      = sys.argv[3] # 2 column file of drug_name - cell_name
    out_file    = sys.argv[4]

    # Extra arguments
    drug_target, sample_mut = None, None
    if len(sys.argv)==6:
    	drug_target = sys.argv[5] # 2 column file of drug_name - gene_target
    elif len(sys.argv)==7:
    	drug_target = sys.argv[5]
    	sample_mut  = sys.argv[6] # matrix of samples x protein_mut

    in_folder   = "/home/ubuntu/"

    output, log = pipeline(in_folder, cell_array, drug_smiles, target, 
    						drug_target=drug_target, sample_mut=sample_mut)
    write_output(out_file, output, log)

    print("Done writing predictions to file\nAppended *.log.pickle file as well")

if __name__ == '__main__':
	run()
