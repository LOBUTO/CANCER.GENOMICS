# ctrp_ridge_cell.py
# Performs Ridge regression on each cell based on gene expression and GSEA-DAE features

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
def small_sample_parse_drug(target_table, drug_table):

	# Parse
    drugs  = list(set(target_table.Compound).intersection(drug_table.index))
    print(len(drugs))

    target_table = target_table.loc[target_table.Compound.isin(drugs)]

    drug_table   = drug_table.loc[drugs]
    target_table = target_table.reset_index(drop=True)

    return target_table, drug_table

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

	# Prediction
	predictions = clf.predict(test_data)

	return predictions

# Load data
ctrp_act  = pd.read_csv("CTRP_FILES/ctrp_act_2.0_true.text", sep="\t")
with open("TEMP_FEATURES/ctrp_drug_feat.pickle", "rb") as handle:
    drug_train, ctrp_dae, drug_dae_log = pickle.load(handle)

# Predict per cell_name
all_cells  = list(set(ctrp_act.cell_name))

main_table = []
count=1.0

for cell in all_cells:
	print(cell)

	train_target  = ctrp_act.query("cell_name==@cell")
	train_feat    = ctrp_dae
	
	train_target, train_feat   = small_sample_parse_drug(train_target, train_feat)

	# We will perform a 5-Fold Cross-validation
	splits      = [list(i) for i in np.array_split(list(train_target.Compound),5)]

	for i in splits:
		train_drugs  = list(set(train_target.Compound).difference(i))

		# Split target data
		split_train  = train_target.query("Compound in @train_drugs")
		split_test   = train_target.query("Compound in @i")
		test_drugs   = split_test.Compound

		split_train_target = np.asarray(split_train.auc)
		split_test_target  = np.asarray(split_test.auc)

		# Split features
		split_train_feat   = train_feat.loc[split_train.Compound]
		split_test_feat    = train_feat.loc[split_test.Compound]

		# Make sure we have enough data
		if split_train_feat.shape[0]>12:

			# Standardize features
			split_train_feat, split_test_feat = scale_standard_multiple(split_train_feat, split_test_feat)

			# Predict
			predictions  = ridge_predict(split_train_feat, split_train_target, split_test_feat)
			drug_pred    = pd.DataFrame({"prediction":predictions, "true":split_test_target,
										"Compound": test_drugs, "cell_name":cell, "feat_type":"smiles_mordred"})

			# Store
			main_table.append(drug_pred)

	# Update count
	count = count + 1
	print("count:", count/(len(all_cells)))

main_table = pd.concat(main_table)

# Store
out_file = "{}ridge_tests_ctrp_drug_smiles_mordred".format(out_folder)
main_table.to_csv(out_file, sep="\t", index=None)
print("Done with ridge testing")

