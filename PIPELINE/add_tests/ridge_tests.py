# Run basic machine learning model: Ridge
# Run on side samples

from sklearn.linear_model import Ridge, RidgeCV
from sklearn.metrics import roc_curve
import cPickle as pickle
import sys
sys.path.insert(0, '/home/ubuntu/GIT/PIPELINE')
from main_pipeline import *
from combat_copy import *

in_folder  = "/home/ubuntu/"
out_folder = "/home/ubuntu/FIGURE_DATA/"

# Minor functions
def small_sample_parse(target_table, exp_table):

	# Parse
    cells  = list(set(target_table.cell_name).intersection(exp_table.index))
    print(len(cells))

    target_table = target_table.loc[target_table.cell_name.isin(cells)]

    exp_table    = exp_table.loc[cells]
    target_table = target_table.reset_index(drop=True)

    return target_table, exp_table

def ridge_predict(train_data, train_target, test_data):

	# Prep modeller
	alpha_ranges = [1e-3, 1e-2, 1e-1, 1, 1e2, 1e3,
					2e3, 2.5e3, 3e3, 3.5e3, 4e3, 
					5e3, 6e3, 6.1e3, 6.15e3, 6.25e3, 6.3e3, 6.4e3, 7e3, 
					7.75e3, 7.9e3, 8e3, 8.1e3, 8.2e3, 8.25e3, 8.3e3, 8.4e3, 8.5e3, 8.75e3, 9e3, 9.25e3, 9.4e3, 9.5e3, 9.6e3, 9.75e3,
					1e4, 1.25e4, 1.4e4, 1.5e4, 1.55e4, 1.58e4, 1.6e4, 1.625e4, 1.65e4, 1.7e4, 1.725e4, 1.74e4, 1.75e4, 1.76e4, 1.78e4, 1.85e4, 
					2e4, 2.25e4, 2.5e4, 3e4, 4e4,  
					0.5e5, 0.75e5, 1e5, 1.25e5, 1.5e5, 
					0.8e6, 0.9e6, 1e6, 1.1e6, 1.2e6, 1.25e6, 1.28e6, 1.3e6, 1.32e6, 1.33e6, 1.34e6, 1.4e6, 1.5e6, 2e6,
					1e7, 1e8, 1e9, 5e9, 1e10, 5e10, 1e11, 1e12, 1e13]
	clf = RidgeCV(alphas=alpha_ranges, 
              normalize=True, cv=None, fit_intercept=False, store_cv_values=True)

	# Fit
	clf.fit(train_data, train_target)
	# print("alpha range:", alpha_ranges)
	# print("CV per alpha:",np.mean(clf.cv_values_, axis=0))
	# print("alpha used:", clf.alpha_)
	# print("fit score:", clf.score(train_data, train_target))

	# Prediction
	predictions = clf.predict(test_data)

	return predictions

# Choice of target variable
choice="pIC50"

# Load gene expresion data for all geeeleher clinical samples (train and test)
doc_exp_train, doc_exp_test       = process_cell_data("GEELEHER/DOCETAXEL/exp_data.txt", in_folder)
bora_exp_train, bora_exp_test     = process_cell_data("GEELEHER/BORTEZOMIB/exp_table_a", in_folder)
borb_exp_train, borb_exp_test     = process_cell_data("GEELEHER/BORTEZOMIB/exp_table_b", in_folder)
cisp_exp_train, cisp_exp_test     = process_cell_data("GEELEHER/CISPLATIN/exp_table", in_folder)

# Load expression feature data for autoencoded gene samples
for i in ["doc", "bora", "borb", "cisp"]:
	for init in [1]:
		for arch in ["2", "2_2", "2_4_8", "4", "4_4", "4_8"]:
			with open("TEMP_FEATURES/{}_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(i, init, arch), "rb") as handle:
				vars()["{}_{}_{}_dae_exp_train".format(i, init, arch)], vars()["{}_{}_{}_dae_exp_test".format(i, init, arch)], vars()["{}_{}_{}_dae_log".format(i, init, arch)] = pickle.load(handle)

for i in ["doc", "bora", "borb", "cisp"]:
	for init in [5]:
		for arch in ["2", "2_2", "2_4_8", "4", "4_4", "4_8"]:
			with open("TEMP_FEATURES/{}_cell_feat_c2setcovertanimoto_testvalid0.4_zscaled_init{}_{}.pickle".format(i, init, arch), "rb") as handle:
				vars()["{}_{}_{}_dae_exp_train".format(i, init, arch)], vars()["{}_{}_{}_dae_exp_test".format(i, init, arch)], vars()["{}_{}_{}_dae_log".format(i, init, arch)] = pickle.load(handle)

# Load target training data
target_table                      = cgp_act_post_process(pd.read_csv("{}CGP_FILES/v17.3_fitted_dose_response.csv".format(in_folder)),
                                                           zscoring=False, choice=choice)

# Predict per compound
main_table = []
drug_dict  = {"doc":"Docetaxel","bora":"Bortezomib","borb":"Bortezomib","cisp":"Cisplatin"}
for init in [1]:
	for arch in ["2", "2_2", "2_4_8", "4", "4_4", "4_8"]:
		for d in ["doc" , "bora", "borb", "cisp"]:
			drug_dict["{}_{}_{}_dae".format(d, init, arch)] = drug_dict[d]

for init in [5]:
	for arch in ["2", "2_2", "2_4_8", "4", "4_4", "4_8"]:
		for d in ["doc", "bora", "borb", "cisp"]:
			drug_dict["{}_{}_{}_dae".format(d, init, arch)] = drug_dict[d]


for ge in drug_dict.keys():
	print(ge)
	drug                = drug_dict[ge]
	dae_error           = None
	print("Drug tested: {}".format(drug))

	# Clean expression and target data
	train_feat			 = cell_gdsc_zhang_rename_index(locals()["{}_exp_train".format(ge)], "gtz")
	ge_table, train_feat = small_sample_parse(target_table, train_feat)

	# Consolidate data
	ge_table    = ge_table.query("Compound==@drug")
	drug_cells  = list(ge_table.cell_name)

	train_feat  = train_feat.loc[drug_cells]
	cell_target = np.asarray(ge_table.value)
	test_feat   = locals()["{}_exp_test".format(ge)]

	# Do we have to separate genes from GSEA-DAE features
	if "dae" in ge:
		dae_error = locals()["{}_log".format(ge)].test_error.values[0] # Update dae_error

		dae_feat  = filter(lambda x: "_gs" in x, list(train_feat.columns))
		genes     = filter(lambda x: "_gs" not in x, list(train_feat.columns))

		train_dae_feat  = np.asarray(train_feat[dae_feat])
		train_gene_feat = np.asarray(train_feat[genes])
		
		test_dae_feat   = np.asarray(test_feat[dae_feat])
		test_gene_feat  = np.asarray(test_feat[genes])

		# Standardize data
		train_dae_feat, test_dae_feat   = scale_standard_multiple(train_dae_feat, test_dae_feat)
		train_gene_feat, test_gene_feat = scale_standard_multiple(train_gene_feat, test_gene_feat)

		# Predict
		dae_predictions  = ridge_predict(train_dae_feat, cell_target, test_dae_feat)
		gene_predictions = ridge_predict(train_gene_feat, cell_target, test_gene_feat)

		# Clean up
		side_predictions = pd.concat([pd.DataFrame({"prediction":dae_predictions,
												   "cell_name":test_feat.index,
												   "kind":"DAE"}),
									  pd.DataFrame({"prediction":gene_predictions,
												   "cell_name":test_feat.index,
												   "kind":"CGC"})])

	train_feat  = np.asarray(train_feat)
	test_feat   = np.asarray(test_feat)

	# Standardize features
	train_feat, test_feat = scale_standard_multiple(train_feat, test_feat)

	# Predict
	predictions  = ridge_predict(train_feat, cell_target, test_feat)
	drug_pred    = pd.DataFrame({"prediction":predictions,
								"cell_name": locals()["{}_exp_test".format(ge)].index,
								"kind":"ALL"})

	if "dae" in ge:
		drug_pred = pd.concat([drug_pred, side_predictions])

	drug_true    = pd.read_csv("/home/ubuntu/GEELEHER/{}/target_table_true.txt".format(drug.upper()),
							  sep="\t")
	drug_pred    = pd.merge(drug_pred, drug_true, on="cell_name")

	k_table = []
	for k in list(set(drug_pred.kind)):

		if choice=="pIC50":
			roc          = roc_auc_score(drug_pred.query("kind==@k")[["target"]], drug_pred.query("kind==@k")[["prediction"]])
			fpr, tpr, th = roc_curve(drug_pred.query("kind==@k")[["target"]], drug_pred.query("kind==@k")[["prediction"]], pos_label=1)
		elif choice=="AUC":
			roc          = roc_auc_score(drug_pred.query("kind==@k")[["target"]], -drug_pred.query("kind==@k")[["prediction"]])
			fpr, tpr, th = roc_curve(drug_pred.query("kind==@k")[["target"]], -drug_pred.query("kind==@k")[["prediction"]], pos_label=1)

		k_table.append(pd.DataFrame({"drug":drug, "source":ge, "kind":k,  "roc":roc, "fpr":fpr, "tpr":tpr, "dae_error":dae_error}))
	k_table = pd.concat(k_table)

	main_table.append(k_table)

main_table = pd.concat(main_table)

# Store
out_file = "{}ridge_tests_{}_test_valid".format(out_folder, choice)
main_table.to_csv(out_file, sep="\t")
print("Done with ridge testing")
