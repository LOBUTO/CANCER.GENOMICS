# main_pipeline.py

import sys
import gc
sys.path.insert(0, '/home/ubuntu/GIT')

from paper_1 import *
from combat_copy import *
from gsea_dae import *

def load_train_exp(in_folder):
	# Loads expression data obtained from GDSC as samplesxgenes
	x = pd.read_csv("{}CGP_FILES/070818_cgp_exp.txt", sep="\t")
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
	train_exp, test_exp = post_process_batch_norm(barch_norm)

	return train_exp.T, test_exp.T

def change_smiles_pubchem_fingerprint(smiles_list):
	# smiles_list: list of smiles strings
	import pubchempy as pcp

	main_list = []
	for i in smiles_list:
		p_id        = pcp.get_cids(identifier=i, namespace="smiles")
		fingerprint = pcp.Compound.from_cid(148124).cactvs_fingerprint
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
	features   = change_smiles_pubchem_fingerprint(all_smiles)

	feat_table = pd.DataFrame(features, index=all_drugs,
							  columns = ["s_%s"%j for j in xrange(len(features))] )

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
	# Extracts canonical smiles from pubchem_ID
	import pubchempy as pcp

	# Load files
	file_in    = "{}pubchem_id.csv".format(in_folder)
	file_in_updated = "{}gdsc_drug_list.csv".format(in_folder)

    x          = pd.read_csv(file_in, header=None, names=["Compound", "count", "ID","N", "DROP"]).drop_duplicates() # SOURCE 1
    x_update   = pd.read_csv(file_in_updated, usecols=[1,5], header=0, names=["Compound", "ID"]).drop_duplicates() # SOURCE 2

    # Clean up
    x["Compound"] = [i.rstrip("'") for i in list(x.Compound)]
    x_update      = x_update.query("ID!='none'").query("ID!='several'")
    x_update      = x_update.dropna()
    
    # Combine with updated
    x_update   = x_update.loc[~x_update.Compound.isin(list(x.Compound))]
    x          = pd.concat([x[["Compound", "ID"]], x_update])

    # Extract smiles using pubchem IDs
    all_ids    = list(x.ID)
	all_drugs  = list(x.Compound)

	smiles     = []
	for i in all_ids:
		smiles.append(pcp.Compound.from_cid(i).canonical_smiles)

	# Clean up
	main_table = pd.DataFrame({"Compound":all_drugs, 
							   "smiles":smiles})
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
	train_fp       = pd.read_csv("{}CGP_FILES/pubchem_smiles_updated.txt".format(in_folder), index_col=0)
	train_fp.index = [Function_drug_name_zhang_to_gdsc(i, "gtz") for i in list(train_fp.index)]
	train_fp       = train_fp.rename_axis("Compound")

	# Extract mordred train data
	train_mordred  = construct_mordred_features(train_smiles)

	# Extract fp and mordred test features
	test_fp        = construct_fingerprint_features(test_smiles)
	test_mordred   = construct_mordred_features(test_smiles)

	# Clean up
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
    exp_genes = exp_genes.intersection(genes)
    
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
	test_var    = pd.read_csv(test_var, sep="\t", names=["Gene", "AA"])
	test_var["effect_name"] = test_var["Gene"] + "_" + gdsc["AA"]
    test_var["value"] = 1
    test_var    = test_var.pivot_table(index="SAMPLE", columns="effect_name", values="value", fill_value=0)

    # Assess compatibility
    common_vars = list(set(train_var.columns).intersection(test_var.columns))

    # Return if common variants found
    if len(common_vars)>2:
    	print("{} common variants found between training and test data".format(len(common_vars)))
    	train_var, test_var = train_var[common_vars], test_var[common_vars]
    	return train_var, test_var
    else:
    	return None, None

def cell_features(cell_exp_train, cell_exp_test, in_folder, CGC=False, test_var=None):
	# Function to extract ALL types of cell-related features (expression, variation)
	# *_exp_train: as panda matrices as samples x genes
	# test_var: 2-column tab-delimited file containing:
	#	- First column:  Gene-name
	#	- Second column: protein mutation as "p.A#B", where:
	#		- A: Original amino acid, B: Final amino acid, # amino acid position
	#		- Can place "*" in place of "B" 

	# Obtain GSEA-DAE features
	train_dae, test_dae, log = gsea_autoencoder(cell_exp_train, cell_exp_test, in_folder,
												gsea="c2setcover", filter=int(10), arch="2")

	# Obtain CGC-GIANT features (if necessary)
	if CGC==True:
		th    = 0.98
		train_cgc, test_cgc = process_giant_features(cell_exp_train, cell_exp_test, in_folder, th)

	# Obtain Variant features (if available)
	if train_var is not None:
		train_var, test_var = process_var_features(test_var, in_folder)

	# Consolidate all features
	train_features = pd.concat([train_dae, train_cgc, train_var], axis=1)
	test_features  = pd.concat([test_dae, test_cgc, test_var], axis=1)

	# Clean up and return
	train_features = train_features.dropna()
	test_features  = test_features.dropna()

	return train_features, test_features

def drug_features(drug_smiles, in_folder):
	 # Function to extract ALL types of drug-related features

	 # Extract smiles-related features first
	 train_fp, train_mordred, test_fp, test_mordred = process_drug_smiles(drug_smiles, in_folder)

	 # Autoencode mordred-based features
	 train_mordred, test_mordred  = feature_dae(train_mordred, )

def pipeline(in_folder, cell_array, drug_smiles, target):
	# cell_array: text array of samples x genes
	# drug_smiles: 2 column file of compounds and smiles, tab-separated
	#	ie. CTRP_FILES/drug_smiles.txt

	# Batch normalize gene expression data - samples x genes
	cell_exp_train, cell_exp_test = process_cell_data(cell_array, in_folder)

	# Extract cell features
	cell_train, cell_test         = cell_features(cell_exp_train, cell_exp_test, in_folder)

	# Extract drug features
	drug_train, drug_test         = drug_features(drug_smiles, in_folder)
	
	# Model

	# Predict

def run():
    cell_array  = sys.argv[1] # cell x genes file, such as ctrp_exp_2.1.txt
    drug_smiles = sys.argv[2] # 2 column file of drug_name - smiles
    target      = sys.argv[3] # 3 column file of drug_name - cell_name - target

    in_folder   = "/home/ubuntu/"
    output      = pipeline(in_folder, cell_array, drug_smiles, target)

if __name__ == '__main__':
	run()
