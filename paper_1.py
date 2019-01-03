# Scripts paper_1.py

import sys
import tensorflow as tf
import pandas as pd
import numpy as np
import random
import pickle
import math
from sklearn.metrics import roc_auc_score
from subprocess import Popen, PIPE
from sklearn import preprocessing
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

# sys.path.insert(0, '/path/to/application/app/folder')

def softmax_matrix(a):
    return(np.exp(a)/np.sum(np.exp(a)))

def process_paper_lobico(lobico_file):
    # lobico_file such as "PAPERS/LORIO/lobico.csv"
    
    # Read in 
    lobico = pd.read_csv(lobico_file, sep=",", na_values=" ")

    lobico = pd.melt(lobico, id_vars="Screened Compounds:")
    lobico.columns = ["cell_name", "Compound", "value"]

    # Clean up
    lobico.value = lobico.value.replace({"S ":"S", "R ":"R"})
    lobico = lobico[lobico.value.notnull()]
    lobico.value = lobico.value.replace({"S":1, "R":0})

    lobico.Compound = [i.rstrip(" ") for i in lobico.Compound]

    return lobico

def process_paper_zhang(feat_file, id_file):
    x = pd.read_csv(feat_file, sep=",")
    print(x.shape)
    y = pd.read_csv(id_file, sep=",", header=None).iloc[:,[0,2]]
    
    x = pd.melt(x, id_vars="Name")
    x.columns = ["ID", "feature", "value"]
    y.columns = ["Compound", "ID"]
    y.Compound = [i.rstrip("'") for i in y.Compound]
    
    x = pd.merge(x, y, on="ID")
    x = x.iloc[:,1:]
    
    return x

def exp_encode_parse(gsea, g_filter, arch, dropkeep, in_folder="/home/ubuntu/"):
    in_folder = "{}GSEA_FILES/RESULTS/".format(in_folder)
    
    main_list = []
    for i in ["train", "valid", "test"]:
        in_file = "{}CGP_{}_{}_deepautoencoder_{}_{}_{}.txt".format(in_folder,gsea,g_filter,arch,dropkeep,i)
        
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

def drug_encode_parse(arch, dropkeep, in_folder="/home/ubuntu/"):
    in_folder = "{}CGP_FILES/RESULTS/".format(in_folder)
    
    main_list = []
    for i in ["train", "valid", "test"]:

        in_file = "{}{}_deepautoencoder_{}_{}_{}.txt".format(in_folder, "cgp_drugs", arch, dropkeep,i)
        main_list.append(pd.read_csv(in_file, sep="\t", index_col=0))
    
    return pd.concat(main_list)

def drug_exp_lobico_parse(lobico, exp, drug):
    
    # Fix for name matching in lobico
    lobico.cell_name = lobico.cell_name.replace({"G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7","Hs633T":"Hs-633T",
                                                 "Hs683":"Hs-683", "Hs766T":"Hs-766T","Hs940-T":"Hs-940-T",
                                                 "HUTU-80":"HuTu-80", "NTERA-2 cl.D1":"NTERA-S-cl-D1",
                                                 "PC-3 [JPC-3]":"PC-3_[JPC-3]", "PE/CA-PJ15":"PE-CA-PJ15",
                                                 "NB-4":"NB4", "G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7"
                                                })
    lobico.Compound  = lobico.Compound.replace({"Zibotentan, ZD4054":"Zibotentan", "PXD101, Belinostat": "Belinostat"}) 
    # Parse
    cells  = list(set(lobico.cell_name).intersection(exp.index))
    drugs  = list(set(lobico.Compound).intersection(drug.index))
    print(len(cells))
    print(len(drugs))
    
    lobico = lobico.loc[lobico.cell_name.isin(cells)]
    lobico = lobico.loc[lobico.Compound.isin(drugs)]
    
    drug   = drug.loc[drugs]
    exp    = exp.loc[cells]
    
    return lobico, exp, drug

def drug_exp_lobico_parse_v2(lobico, exp, mut, drug, replacement="ztg", binary_target=False):
    # Replacement of compound names for compatibilitiy as: zhang_to_gdsc(ztg) or gdsc_to_zhang(gtz)

    # Fix for name matching in lobico
    lobico.cell_name = lobico.cell_name.replace({"G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7","Hs633T":"Hs-633T",
                                                 "Hs683":"Hs-683", "Hs766T":"Hs-766T","Hs940-T":"Hs-940-T",
                                                 "HUTU-80":"HuTu-80", "NTERA-2 cl.D1":"NTERA-S-cl-D1",
                                                 "PC-3 [JPC-3]":"PC-3_[JPC-3]", "PE/CA-PJ15":"PE-CA-PJ15",
                                                 "NB-4":"NB4", "G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7"
                                                })
    lobico.Compound  = lobico.Compound.replace({"Zibotentan, ZD4054":"Zibotentan", "PXD101, Belinostat": "Belinostat"})

    # Added - Replacement added due to published compound features matrix drug names on Zhang et al. 2018
    gdsc_to_zhang    = {"AICA Ribonucleotide":"AICAR","BAY-61-3606":"BAY 61-3606", "BX795":"BX-795",
    					"CCT-018159":"CCT018159", "EHT-1864": "EHT 1864", "GSK1904529A": "GSK-1904529A",
    					"GW441756":"GW 441756", "HG6-64-1":"HG-6-64-1", "Nutlin-3a (-)":"Nutlin-3a",
    					"PD0325901":"PD-0325901", "PD173074":"PD-173074", "PLX-4720":"PLX4720",
    					"SB216763":"SB 216763", "SB505124":"SB-505124", "SL0101":"SL 0101-1",
    					"XAV939":"XAV 939", "ZM447439":"ZM-447439",
    					"Tanespimycin":"17-AAG", "Wee1 Inhibitor":"681640", "Navitoclax":"ABT-263",
    					"Linifanib":"ABT-869", "Veliparib":"ABT-888", "Quizartinib":"AC220",
    					"Rucaparib":"AG-014699", "Motesanib":"AMG-706", "Ponatinib":"AP-24534",
    					"Tretinoin":"ATRA", "Luminespib":"AUY922", "Tivozanib":"AV-951", "Saracatinib":"AZD-0530",
    					"Selumetinib":"AZD6244", "Talazoparib":"BMN-673", "Avagacestat":"BMS-708163",
    					"Idelalisib":"CAL-101", "Lestaurtinib":"CEP-701", "Alectinib":"CH5424802",
    					"Pelitinib":"EKB-569", "Selisistat":"EX-527",
    					"Daporinad":"FK866", "Pictilisib":"GDC0941", "Omipalisib":"GSK2126458",
    					"Serdemetan":"JNJ-26854165", "Dacinostat":"LAQ824", "Enzastaurin":"LY317615",
    					"Pevonedistat":"MLN4924", "Amuvatinib":"MP470", "Entinostat":"MS-275",
    					"Dactolisib":"NVP-BEZ235", "Linsitinib":"OSI-906", "Palbociclib":"PD-0332991",
    					"Refametinib":"RDEA119", "Seliciclib":"Roscovitine", "Ispinesib Mesylate":"SB-715992",
    					"Fedratinib":"TG101348", "Tozasertib":"VX-680", "Cabozantinib":"XL-184",
    					"Foretinib":"XL-880", "Sepantronium bromide":"YM155"}

    zhang_to_gdsc    = {y:x for x,y in gdsc_to_zhang.iteritems()}

    if replacement=="gtz":
    	print("renaming GDSC compounds names to published Zhang 2018 names")
    	lobico.Compound  = lobico.Compound.replace(gdsc_to_zhang)
    elif replacement=="ztg":
    	print("renaming Zhang published names to original GDSC names")
    	drug = drug.rename(index=zhang_to_gdsc)

    # Parse
    cells  = set(lobico.cell_name).intersection(exp.index)
    cells  = list(set(cells).intersection(mut.index))
    drugs  = list(set(lobico.Compound).intersection(drug.index))
    print(len(cells))
    print(len(drugs))

    lobico = lobico.loc[lobico.cell_name.isin(cells)]
    lobico = lobico.loc[lobico.Compound.isin(drugs)]

    drug   = drug.loc[drugs]
    exp    = exp.loc[cells]
    mut    = mut.loc[cells]

    return lobico, exp, mut, drug

def cell_gdsc_zhang_rename_index(cell_data, replacement="ztg"):
	# Renames index of matrix based on zhang or gdsc nomenclature

	gdsc_to_zhang = {"G-292 Clone A141B1":"G-292-Clone-A141B1",
					 "Hep 3B2_1-7":"Hep3B2-1-7","Hs633T":"Hs-633T",
					 "Hs683":"Hs-683", "Hs766T":"Hs-766T","Hs940-T":"Hs-940-T",
					 "HUTU-80":"HuTu-80", "NTERA-2 cl.D1":"NTERA-S-cl-D1",
					 "PC-3 [JPC-3]":"PC-3_[JPC-3]", "PE/CA-PJ15":"PE-CA-PJ15",
					 "NB-4":"NB4", "G-292 Clone A141B1":"G-292-Clone-A141B1",
					 "Hep 3B2_1-7":"Hep3B2-1-7"}

	zhang_to_gdsc = {y:x for x,y in gdsc_to_zhang.iteritems()}

	if replacement=="gtz":
		print("renaming GDSC compounds names to published Zhang 2018 names")
		cell_data = cell_data.rename(index=gdsc_to_zhang)
	elif replacement=="ztg":
		print("renaming Zhang published names to original GDSC names")
		cell_data = cell_data.rename(index=zhang_to_gdsc)

	return cell_data

def drug_gdsc_zhang_rename_index(drug_data, replacement="ztg"):
	# Renames index of matrix based on zhang or gdsc nomenclature

	# Clean up first
	drug_data = drug_data.rename(index={"Zibotentan, ZD4054":"Zibotentan", "PXD101, Belinostat": "Belinostat"})

	# Build renaming dict
	gdsc_to_zhang    = {"AICA Ribonucleotide":"AICAR","BAY-61-3606":"BAY 61-3606", "BX795":"BX-795",
						"CCT-018159":"CCT018159", "EHT-1864": "EHT 1864", "GSK1904529A": "GSK-1904529A",
						"GW441756":"GW 441756", "HG6-64-1":"HG-6-64-1", "Nutlin-3a (-)":"Nutlin-3a",
						"PD0325901":"PD-0325901", "PD173074":"PD-173074", "PLX-4720":"PLX4720",
						"SB216763":"SB 216763", "SB505124":"SB-505124", "SL0101":"SL 0101-1",
						"XAV939":"XAV 939", "ZM447439":"ZM-447439",
						"Tanespimycin":"17-AAG", "Wee1 Inhibitor":"681640", "Navitoclax":"ABT-263",
						"Linifanib":"ABT-869", "Veliparib":"ABT-888", "Quizartinib":"AC220",
						"Rucaparib":"AG-014699", "Motesanib":"AMG-706", "Ponatinib":"AP-24534",
						"Tretinoin":"ATRA", "Luminespib":"AUY922", "Tivozanib":"AV-951", "Saracatinib":"AZD-0530",
						"Selumetinib":"AZD6244", "Talazoparib":"BMN-673", "Avagacestat":"BMS-708163",
						"Idelalisib":"CAL-101", "Lestaurtinib":"CEP-701", "Alectinib":"CH5424802",
						"Pelitinib":"EKB-569", "Selisistat":"EX-527",
						"Daporinad":"FK866", "Pictilisib":"GDC0941", "Omipalisib":"GSK2126458",
						"Serdemetan":"JNJ-26854165", "Dacinostat":"LAQ824", "Enzastaurin":"LY317615",
						"Pevonedistat":"MLN4924", "Amuvatinib":"MP470", "Entinostat":"MS-275",
						"Dactolisib":"NVP-BEZ235", "Linsitinib":"OSI-906", "Palbociclib":"PD-0332991",
						"Refametinib":"RDEA119", "Seliciclib":"Roscovitine", "Ispinesib Mesylate":"SB-715992",
						"Fedratinib":"TG101348", "Tozasertib":"VX-680", "Cabozantinib":"XL-184",
						"Foretinib":"XL-880", "Sepantronium bromide":"YM155"}

	zhang_to_gdsc    = {y:x for x,y in gdsc_to_zhang.iteritems()}

	if replacement=="gtz":
		print("renaming GDSC compounds names to published Zhang 2018 names")
		drug_data = drug_data.rename(index=gdsc_to_zhang)
	elif replacement=="ztg":
		print("renaming Zhang published names to original GDSC names")
		drug_data = drug_data.rename(index=zhang_to_gdsc)

	return drug_data

def drug_exp_lobico_parse_v3(lobico, exp, drug, replacement="ztg", binary_target=False):
    # Replacement of compound names for compatibilitiy as: zhang_to_gdsc(ztg) or gdsc_to_zhang(gtz)
    # Updated so that it is fed only two sources of information: expression and drug

    # Fix for name matching in lobico
    lobico.cell_name = lobico.cell_name.replace({"G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7","Hs633T":"Hs-633T",
                                                 "Hs683":"Hs-683", "Hs766T":"Hs-766T","Hs940-T":"Hs-940-T",
                                                 "HUTU-80":"HuTu-80", "NTERA-2 cl.D1":"NTERA-S-cl-D1",
                                                 "PC-3 [JPC-3]":"PC-3_[JPC-3]", "PE/CA-PJ15":"PE-CA-PJ15",
                                                 "NB-4":"NB4", "G-292 Clone A141B1":"G-292-Clone-A141B1",
                                                 "Hep 3B2_1-7":"Hep3B2-1-7"
                                                })
    lobico.Compound  = lobico.Compound.replace({"Zibotentan, ZD4054":"Zibotentan", "PXD101, Belinostat": "Belinostat"})

    # Added - Replacement added due to published compound features matrix drug names on Zhang et al. 2018
    gdsc_to_zhang    = {"AICA Ribonucleotide":"AICAR","BAY-61-3606":"BAY 61-3606", "BX795":"BX-795",
    					"CCT-018159":"CCT018159", "EHT-1864": "EHT 1864", "GSK1904529A": "GSK-1904529A",
    					"GW441756":"GW 441756", "HG6-64-1":"HG-6-64-1", "Nutlin-3a (-)":"Nutlin-3a",
    					"PD0325901":"PD-0325901", "PD173074":"PD-173074", "PLX-4720":"PLX4720",
    					"SB216763":"SB 216763", "SB505124":"SB-505124", "SL0101":"SL 0101-1",
    					"XAV939":"XAV 939", "ZM447439":"ZM-447439",
    					"Tanespimycin":"17-AAG", "Wee1 Inhibitor":"681640", "Navitoclax":"ABT-263",
    					"Linifanib":"ABT-869", "Veliparib":"ABT-888", "Quizartinib":"AC220",
    					"Rucaparib":"AG-014699", "Motesanib":"AMG-706", "Ponatinib":"AP-24534",
    					"Tretinoin":"ATRA", "Luminespib":"AUY922", "Tivozanib":"AV-951", "Saracatinib":"AZD-0530",
    					"Selumetinib":"AZD6244", "Talazoparib":"BMN-673", "Avagacestat":"BMS-708163",
    					"Idelalisib":"CAL-101", "Lestaurtinib":"CEP-701", "Alectinib":"CH5424802",
    					"Pelitinib":"EKB-569", "Selisistat":"EX-527",
    					"Daporinad":"FK866", "Pictilisib":"GDC0941", "Omipalisib":"GSK2126458",
    					"Serdemetan":"JNJ-26854165", "Dacinostat":"LAQ824", "Enzastaurin":"LY317615",
    					"Pevonedistat":"MLN4924", "Amuvatinib":"MP470", "Entinostat":"MS-275",
    					"Dactolisib":"NVP-BEZ235", "Linsitinib":"OSI-906", "Palbociclib":"PD-0332991",
    					"Refametinib":"RDEA119", "Seliciclib":"Roscovitine", "Ispinesib Mesylate":"SB-715992",
    					"Fedratinib":"TG101348", "Tozasertib":"VX-680", "Cabozantinib":"XL-184",
    					"Foretinib":"XL-880", "Sepantronium bromide":"YM155"}

    zhang_to_gdsc    = {y:x for x,y in gdsc_to_zhang.iteritems()}

    if replacement=="gtz":
    	print("renaming GDSC compounds names to published Zhang 2018 names")
    	lobico.Compound  = lobico.Compound.replace(gdsc_to_zhang)
    elif replacement=="ztg":
    	print("renaming Zhang published names to original GDSC names")
    	drug = drug.rename(index=zhang_to_gdsc)

    # Parse
    cells  = list(set(lobico.cell_name).intersection(exp.index))
    drugs  = list(set(lobico.Compound).intersection(drug.index))
    print(len(cells))
    print(len(drugs))

    lobico = lobico.loc[lobico.cell_name.isin(cells)]
    lobico = lobico.loc[lobico.Compound.isin(drugs)]

    drug   = drug.loc[drugs]
    exp    = exp.loc[cells]
    lobico = lobico.reset_index(drop=True)

    return lobico, exp, drug

def cleans_ctrp_act(ctrp_act_file):
	print("Cleaning up multiple entries in CTRP")
	# Fix for compounds that have multiple measurements
	# Cleans up ctrp

	ctrp_act      = pd.read_csv(ctrp_act_file, sep="\t")
	ctrp_act["mean_auc"] = ctrp_act.groupby(["Compound","cell_name"])["auc"].transform(np.mean)

	ctrp_act["auc"] = ctrp_act["mean_auc"]
	ctrp_act      = ctrp_act.drop("mean_auc", axis=1)
	ctrp_act["neg_auc"] = -ctrp_act["auc"]
	ctrp_act      = ctrp_act.drop_duplicates()

	return(ctrp_act)

def process_cgp_multiple_conc(cgp_act):
	print("Cleaning up multiple entries in GDSC")
	# Fix for compounds that have multiple measurements
	x             = cgp_act.groupby(["DRUG_NAME", "DRUG_ID"])["DRUG_NAME"].count().reset_index(name="COUNT")
	m_drugs       = list((x.DRUG_NAME.value_counts().reset_index(name="cell_count")
						.query("cell_count>1"))["index"])

	x_filter     = cgp_act.loc[~cgp_act.DRUG_NAME.isin(m_drugs)]

	x_fix        = x.loc[x.DRUG_NAME.isin(m_drugs)]
	x_fix["mak"] = x_fix.groupby("DRUG_NAME")["COUNT"].transform(max)
	x_fix        = x_fix.query("COUNT==mak")

	# Clean up
	x_fix        = pd.merge(cgp_act[["CELL_LINE_NAME", "DRUG_NAME","LN_IC50", "AUC", "DRUG_ID"]],
							x_fix[["DRUG_NAME", "DRUG_ID"]])

	x            = pd.concat([x_filter[["DRUG_NAME", "CELL_LINE_NAME", "LN_IC50", "AUC"]],
								x_fix[["DRUG_NAME", "CELL_LINE_NAME", "LN_IC50", "AUC"]]])
	return(x)

def cgp_act_post_process(cgp_act, zscoring=False, choice="pIC50", rebalance=False):
	# Calculate the pIC50 value across all drugs
	from scipy.stats import zscore

	# Load file
	cgp_act          = process_cgp_multiple_conc(cgp_act)

	# Clean up
	cgp_act["DRUG_NAME"] = [i.rstrip(" ") for i in list(cgp_act.DRUG_NAME)]

	if choice=="pIC50":
		cgp_act["value"] = -np.log(np.exp(cgp_act.LN_IC50))
	elif choice=="AUC":
		cgp_act["value"] = cgp_act.AUC

	cgp_act          = cgp_act[["CELL_LINE_NAME", "DRUG_NAME", "value"]]
	cgp_act.columns  = ["cell_name", "Compound", "value"]

	# Do we need to rebalance?
	if rebalance==True:
		print("Rebalancing data set")
		max_number_cells = float(list(cgp_act.Compound.value_counts().reset_index(name="n").sort_values("n", ascending=False)["n"])[0])

		all_compounds    = list(set(cgp_act.Compound))
		balanced_cgp     = []
		for c in all_compounds:
			
			query_c = cgp_act.query("Compound==@c")
			if query_c.shape[0] < max_number_cells:
				
				balanced_cgp.append(query_c.sample(frac=max_number_cells/query_c.shape[0], replace=True, random_state=1))

		cgp_act = pd.concat(balanced_cgp).reset_index(drop=True)

	# Do we need to z-scale the target variable?
	if zscoring==True:
		print("Z-scoring target value")
		cgp_act_z = [pd.DataFrame({"cell_name":j.cell_name, "Compound":i, "value":zscore(j.value)}) for
					i,j in cgp_act.groupby("Compound")]
		cgp_act_z = pd.concat(cgp_act_z)
		return(cgp_act_z)

	else:
		return(cgp_act)

# DEEP LEARNING TENSORFLOW FUNCTIONS
def weight_variable(shape):
	initial = tf.truncated_normal(shape, stddev=0.1) #Normal distributionwith std 0.1
	return tf.Variable(initial)

def bias_variable(shape):
	initial = tf.constant(0.1, shape=shape) #Uniform distribution of 0.1
	return tf.Variable(initial)

def lrelu(x,alpha=0.1):
	return tf.maximum(alpha*x,x)

def valid_indexes(line_count):

    random.seed(1234) #NOTE Remove for total random (rnd)
    temp_index  = random.sample(xrange(line_count), line_count) #Randomize at first
    valid_index = list(np.array_split(temp_index, 5)[4]) #Validate on 10%

    train_index = list(set(temp_index) - set(valid_index)) #Train on 90%
    print("Train samples: %s, Valid samples: %s"%(len(train_index), len(valid_index)) )

    return train_index, valid_index

def data_split(samples_list):
	random.seed(1234)
	temp_index  = random.sample(samples_list, len(samples_list))
	valid_index = list(np.array_split(temp_index, 10)[9])
	train_index = list(set(temp_index) - set(valid_index)) #Train on 80%

	print("Train samples: %s, Valid samples: %s"%(len(train_index), len(valid_index)) )
	return train_index, valid_index

def drug_feat_target_cv_split(target_table, feat_table, drug, cells):
	# Cross-validation split based on set of cells (ie. 10%) per drug

    # Index based on training sample
    all_index    = list(target_table.index)
    valid_index  = list(target_table.loc[(target_table.Compound==drug) & (target_table.cell_name.isin(cells))].index)
    print("valid_index: ", valid_index)
    train_index  = list(set(all_index) - set(valid_index))

    train_feat   = np.asarray(feat_table.loc[train_index])
    valid_feat   = np.asarray(feat_table.loc[valid_index])
    
    train_target = np.asarray(target_table.loc[train_index].value)
    valid_target = np.asarray(target_table.loc[valid_index].value)

    # Return true_values
    true_values  = target_table.loc[valid_index]["value"].values

    return train_feat, valid_feat, train_target, valid_target, true_values

def drug_feat_target_split(target_table, feat_table, drug):
	# Split leaving target drug along with all its target cells  on valid set

	all_index    = list(target_table.index)
	valid_index  = list(target_table.loc[target_table.Compound==drug].index)
	train_index  = list(set(all_index) - set(valid_index))

	train_feat   = np.asarray(feat_table.loc[train_index])
	valid_feat   = np.asarray(feat_table.loc[valid_index])

	train_target = np.asarray(target_table.loc[train_index].value)
	valid_target = np.asarray(target_table.loc[valid_index].value)

	# Return true_values
	true_values  = target_table.loc[valid_index]["value"].values

	return train_feat, valid_feat, train_target, valid_target, true_values

def data_split_drug(target_table, drug):
	# Same purpose as drug_feat_target_split() but it produces only indices rather than data
	# DEPRECATED 
	target_table = target_table.loc[target_table.Compound!=drug]

	drugs        = set(target_table.Compound)
	random.seed(1234)
	train_drugs  = random.sample(drugs, int(math.ceil(len(drugs)*0.8)) )
	valid_drugs  = list(drugs - set(train_drugs))

	train_index  = target_table.loc[target_table.Compound.isin(train_drugs)].index
	valid_index  = target_table.loc[target_table.Compound.isin(valid_drugs)].index

	random.seed(1234)
	train_index  = random.sample(train_index, len(train_index))
	random.seed(1234)
	valid_index  = random.sample(valid_index, len(valid_index))

	return train_index, valid_index

def data_split_drug_self(target_table, drug):
	# Similar to data_split_drug(), but we resample valid from train (subsetting)
	# Applied whe we actually want to train with all of the data

	target_table = target_table.loc[target_table.Compound!=drug]

	drugs        = set(target_table.Compound)
	random.seed(1234)
	train_drugs  = drugs #TRAINING WITH ALL
	valid_drugs  = random.sample(drugs, int(math.ceil(len(drugs)*0.1))) #VALID WITH 10% OF THEM (SUBSET)

	train_index  = target_table.loc[target_table.Compound.isin(train_drugs)].index
	valid_index  = target_table.loc[target_table.Compound.isin(valid_drugs)].index

	random.seed(1234)
	train_index  = random.sample(train_index, len(train_index))
	random.seed(1234)
	valid_index  = random.sample(valid_index, len(valid_index))

	return train_index, valid_index

def train_valid_test_split(line_count):
	# Splits 80/10/10
	random.seed(1234) #NOTE Remove for total random (rnd)

	temp_index  = random.sample(xrange(line_count), line_count) #Randomize at first
	temp_index  = np.array_split(temp_index, 5)

	validation_index        = np.array_split(temp_index[4], 2)
	valid_index, test_index = list(validation_index[0]), list(validation_index[1]) #Validate on 10% and test on 10%
	train_index = list(np.concatenate(temp_index[:4])) #Train on 80%

	return train_index, valid_index, test_index

def scale_0_1_multiple(x, y, z=None):
	# Will scale between 0-1 by column
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per gene across dataset 0-1: multiple")
	min_max_scaler = preprocessing.MinMaxScaler()
	x = min_max_scaler.fit_transform(x)

	if (y is not None):
		y = min_max_scaler.transform(y)

	if (z is not None):
		z = min_max_scaler.transform(z)
		return x, y, z
	else:
		return x, y

def scale_standard_multiple(x, y, z=None):
	# Will zero-center and divide by standard dev.
	# Scaling parameters will be obtained for train(x), and then applied to train and test data(y)

	print("Scaling per gene across dataset standard: multiple")
	scaler = preprocessing.StandardScaler()
	x = scaler.fit_transform(x)

	if (y is not None):
		y = scaler.transform(y)

	if (z is not None):
		z = scaler.transform(z)
		return x, y, z
	else:
		return x, y

def arch_layers(start, arch):
	arch   = arch.split("_")

	layers = [start]
	print(layers)
	for i in arch:
		i = float(i)
		i = start/i
		layers.append(int(i))

	return layers

def dense_batch_lrelu(x, phase, scope):
	# x: matmul(a,b), no bias needed
	with tf.variable_scope(scope):
		h = tf.contrib.layers.batch_norm(x,
			center=True, scale=True,
			is_training=phase,
			scope='bn')

		return lrelu(h, alpha=0.2)

def dense_lrelu_batch(x, phase, scope):
	# x: matmul(a,b), no bias needed
	with tf.variable_scope(scope):
		l = lrelu(x, alpha=0.2)

		h = tf.contrib.layers.batch_norm(l,
			center=True, scale=True,
			is_training=phase,
			scope='bn')
		return h

def binary_one_hot(bin_target):
    # print(bin_target)
    # takes in a numpy array of dim (n,)
    # converts binary vector to one hot
    n_samples = bin_target.shape[0]

    b = np.zeros((n_samples, 2))
    b[np.arange(n_samples), bin_target] = 1

    # print(b)
    return(b)

def entrez_to_hugo(entrez_list, email):
	# Converts entrez IDs to HUGO symbols
    from Bio import Entrez
    Entrez.email = email

    # Split into chunks (max for Bio is 9999)
    a       = [str(i) for i in list(entrez_list)]
    a_parts = [a[i:(i+5000)] for i in xrange(0, len(a), 5000)]
    
    annot_list = []
    for p in a_parts:
        request  = Entrez.epost("gene",id=",".join(p))
        result   = Entrez.read(request)
        
        webEnv   = result["WebEnv"]
        queryKey = result["QueryKey"]
        data     = Entrez.esummary(db="gene", webenv=webEnv, query_key =queryKey)
        annotations = Entrez.read(data)
        
        for gene_data in annotations["DocumentSummarySet"]["DocumentSummary"]:
            entrez_id   = str(gene_data.attributes["uid"])
            gene_symbol = gene_data["NomenclatureSymbol"]
            annot_list.append([entrez_id, gene_symbol])
    
    return(pd.DataFrame(annot_list, columns=["entrez", "gene_name"]))

def process_drug_target(file_in):
	# Process drug target file from GSDC (Screened_Compounds.csv) so that correct HUGO symbol targets are found for all compounds
	drug_target = pd.read_csv(file_in, sep=",")
	drug_target = {drug_target["DRUG_NAME"][i].rstrip(" "):drug_target["TARGET"][i].split(",") for i in xrange(drug_target.shape[0])}
	drug_target = pd.DataFrame([[j,i] for j in drug_target.keys() for i in drug_target[j]], columns=["Compound", "target"])

	drug_target["target"] = [i.lstrip().rstrip() for i in drug_target["target"]]

	# Fixes for common gene names
	drug_target.target = drug_target.target.replace({"ABL":"ABL1", "ABL(T315I)":"ABL1",
													 "ALK4":"ACVR1B", "ALK5":"TGFBR1", "AMPK agonist":"PRKAA1",
													 "Amyloid beta20":"APP", "Amyloid beta40":"APP",
													 "Angiopoietin-1 receptor":"ANGPT1", "Anthracycline":"TOP2B",
													 "BCL-W":"BCL2L2", "BCL-XL":"BCL2L1", "CAMK2":"CAMK2A",
													 "DNAPK":"PRKDC", "Dihydrofolate reductase (DHFR)":"DHFR",
													 "ERK1":"MAPK3", "ERK2":"MAPK1","ERK5":"MAPK7",
													 "Endothelin-1 receptor (EDNRA)":"EDNRA",
													 "FAK":"PTK2","FAK2":"PTK2B","FXR":"NR1H4",
													 "Farnesyl-transferase (FNTA)":"FNTA",
													 "RAR":"RARA", "HIF-PH":"EGLN1",
													 "IKK1":"CHUK", "IKK2":"IKBKB", "IKKB":"IKBKB",
													 "IR":"INSR",
													 "JNK":"MAPK8", "JNK1":"MAPK8", "JNK2":"MAPK9","JNK3":"MAPK10",
													 "LOK":"STK10", "MEK1":"MAP2K1", "MEK2":"MAP2K2", "MEK5":"MAP2K5",
													 "MPS1":"TTK",
													 "PI3K (Class 1)":"PIK3CA", "PI3K (class 1)":"PIK3CA", "PI3Kalpha":"PIK3CA",
													 "PB1":"PBRM1", "PDK1 (PDPK1)":"PDPK1",
													 "PI3Kbeta":"PIK3CB", "PI3Kdelta":"PIK3CD", "PI3Kgamma":"PIK3CG",
													 "PKC":"PRKCA", "PKCB":"PRKCB",
													 "PPARdelta":"PPARD", "PPARgamma":"PPARG", "PYKFYVE":"PIKFYVE",
													 "Procaspase-3":"CASP3", "Procaspase-7":"CASP7",
													 "RNA helicase A":"DHX9",
													 "RON":"MST1R", "RSK":"RPS6KA1", "RSK2":"RPS6KA3",
													 "S6K1":"RPS6KB1", "SHP-1 (PTPN6)":"PTPN6", "SHP-2 (PTPN11)":"PTPN11",
													 "SPRK1":"SRPK1", "TAK":"MAP3K7", "TAK1":"MAP3K7", "TIE2":"TEK",
													 "TNKS1":"TNKS", "TOP2":"TOP2A",
													 "VEGFR":"FLT1", "VEGFR1":"FLT1", "VEGFR2":"KDR", "VEGFR3":"FLT4", "VEGFR3/FLT4":"FLT4",
													 "c-FGR":"FGR", "gamma-secretase":"APH1A",
													 "p38":"MAPK14", "p38alpha":"MAPK14", "p38beta":"MAPK11"})

	# Adds specific targets based on KEGG
	new_set     = {"Cytarabine":["POLA1", "POLA2", "POLB", "POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POL3", "POLE4", "POLG", "POLG2", "POLI", "POLK", "POLL"],
				   "Methotrexate":["DHFR", "DHFRP1", "DYR", "DHFR2", "DHFRL1", "DHFRP4"],
				   "Gemcitabine":["RRM1", "R1", "RIR1", "RR1"],
				   "5-Fluorouracil":["TYMS", "HST422", "TMS", "TS"],
				   "GNF-2":["BCR", "ABL1"],
				   "XMD8-85":["BRD2", "BRD3", "BRD4"],
				   "Phenformin":["PRKAA1", "PRKAA2"],
				   "Temozolomide":["DNA"],
				   "Mitomycin-C":["DNA"],
				   "Dasatinib":["EFNA1","EFNA2", "EFNA3","EFNA4","EFNA5","EFNB1","EFNB2","EFNB3"],
				   "UNC0638":["EHMT2","EHMT1"],
				   "VNLG/124":["HDAC1","HDAC2","HDAC3","HDAC8"],
				   "Vorinostat":["HDAC1","HDAC2","HDAC3","HDAC8","HDAC4","HDAC5","HDAC7","HDAC9","HDAC6","HDAC10","HDAC11"],
				   "CUDC-101":["HDAC1","HDAC2","HDAC3","HDAC8","HDAC4","HDAC5","HDAC7","HDAC9","HDAC6","HDAC10","HDAC11", "HD1", "HD2"],
				   "Tanespimycin":["HSP90B1","HSP90AB1","HSP90AA1", "EL52", "HEL-S-65p", "HSP86", "HSP89A", "HSP90A", "HSP90N", 
				   					"HSPC1", "HSPCA", "HSPCAL1", "HSPCAL4", "HSPN", "LAP2", "D6S182", "HSP84","HSPC2","HSPCB",
				   					"ECGP", "GP96", "GRP94", "HEL35", "TRA1"],
				   	"Elesclomol":["HSP90B1","HSP90AB1","HSP90AA1", "HSPA2"],
				   	"CCT-018159":["HSP90B1","HSP90AB1","HSP90AA1"],
				   	"Luminespib":["HSP90B1","HSP90AB1","HSP90AA1", "EL52", "HEL-S-65p", "HSP86", "HSP89A", "HSP90A", "HSP90N", 
				   					"HSPC1", "HSPCA", "HSPCAL1", "HSPCAL4", "HSPN", "LAP2", "D6S182", "HSP84","HSPC2","HSPCB",
				   					"ECGP", "GP96", "GRP94", "HEL35", "TRA1"],
				   	"SNX-2112":["HSP90B1","HSP90AB1","HSP90AA1"],
				   	"BX795":["CHUK", "IKBKB", "IKBKG","PDK1"],
				   	"Piperlongumine":["TXNRD1", "XPO1", "CYP2D6", "AKT1", "MTOR"],
				   	"Ispinesib Mesylate":["KIF11","EG5", "HKSP", "KNSL1","MCLMR", "TRIP5"],
				   	"T0901317": ["NR1H3","NR1H4","ABCA1","CD36","RORA","RORC","NR1I2"],
				   	"Rapamycin":["MTOR", "AKT1S1", "MLST8", "RPTOR"],
				   	"Omipalisib":["MTOR", "AKT1S1", "MLST8", "RPTOR", "MAPKAP1", "MLST8", "PRR5", "RICTOR"],
				   	"OSI-027":["MTOR", "AKT1S1", "MLST8", "RPTOR", "MAPKAP1", "MLST8", "PRR5", "RICTOR"],
				   	"Dactolisib":["MTOR", "AKT1S1", "MLST8", "RPTOR", "MAPKAP1", "MLST8", "PRR5", "RICTOR"],
				   	"AZD8055":["MTOR", "AKT1S1", "MLST8", "RPTOR", "MAPKAP1", "MLST8", "PRR5", "RICTOR"],
				   	"Vinorelbine":["TUBB3", "CFEOM3", "CFEOM3A", "FEOM3", "TUBB4", "TUBB4A", "DYT4",
				   				   "TUBB4B","LCAEOD","TUBB2","TUBB2C", "TUBB", "CSCSC1", "M40", "TUBB1", "TUBB5",
				   				   "TUBB8", "OOMD", "OOMD2", "TUBB2B", "PMGYSA", "TUBB2A", "TUBB1",
				   				   "TUBB6", "FPVEPD"],
				   	"Vinblastine":["TUBB3", "CFEOM3", "CFEOM3A", "FEOM3", "TUBB4", "TUBB4A", "DYT4",
				   				   "TUBB4B","LCAEOD","TUBB2","TUBB2C", "TUBB", "CSCSC1", "M40", "TUBB1", "TUBB5",
				   				   "TUBB8", "OOMD", "OOMD2", "TUBB2B", "PMGYSA", "TUBB2A", "TUBB1",
				   				   "TUBB6", "FPVEPD"],
				   	"ICL1100013":["NMT1", "NMT2"],
				   	"Pevonedistat":["NAE1", "APPBP1", "HPP1"],
				   	"Sunitinib":["CSF1", "MCSF", "PDGFRA"],
				   	"Midostaurin":["FLT3", "FLK2", "STK1"],
				   	"MG-132":["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMA8", 
				   			  "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10"],
				   	"Bortezomib":["PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMA8",
				   				  "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10"],
				   	"Sorafenib":["RAF1", "CMD1NN", "CRAF", "NS5"],
				   	"CX-5461":["POLR1A", "POLR1B", "POLR1C", "POLR1D", "POLR1E"],
				   	"Y-39983":["ROCK1", "ROCK2"],
				   	"Bexarotene":["RXRA", "NR2B1", "RXRB", "DAUDI6", "H-2RIIBP", "NR2B2", "RXRG", "NR2B3", "RXRC"],
				   	"Tretinoin":["PML", "MYL", "PP8675", "RNF71", "TRIM19", "RARA", "NR1B1", "RAR"],
				   	"Thapsigargin":["ATP2A1", "ATP2A2", "ATP2A3"],
				   	"rTRAIL":["TNFRSF10A", "TNFRSF10B"],
				   	"Shikonin":["PKM", "ANO1"]}
	new_set    = pd.DataFrame([[j,i] for j in new_set.keys() for i in new_set[j]], columns=["Compound", "target"])


	# Make sure that all targets of BX795 are complete
	drug_target = pd.concat([drug_target, new_set])

	return(drug_target)

def string_to_hugo(string_list, alias_file):
	# alias_file like "STRING/9606.protein.aliases.v10.5.txt"

	# Process hugo alias
	alias = pd.read_csv(alias_file, sep="\t", skiprows=1)

	alias.columns = ["ensembl", "alias", "source"]
	alias = alias.loc[alias.source.apply(lambda x: (" Ensembl_HGNC " in x) | ("BioMart_HUGO" in x) )]
	alias = alias[["ensembl", "alias"]]

	# Find names
	alias = alias.loc[alias.ensembl.isin(string_list)]

	# Return
	return alias

def string_pp_to_hugo(pp_file, alias_file):
    # pp_file like "/home/ubuntu/STRING/9606.protein.links.v10.5.txt"
    # alias_file like "STRING/9606.protein.aliases.v10.5.txt"
    
    pp    = pd.read_csv(pp_file, sep=" ")
    alias = pd.read_csv(alias_file, sep="\t", skiprows=1)
    
    alias.columns = ["ensembl", "alias", "source"]
    alias = alias.loc[alias.source.apply(lambda x: (" Ensembl_HGNC " in x) | ("BioMart_HUGO" in x) )]
    alias = alias[["ensembl", "alias"]]
    
    pp    = pd.merge(pp, alias, left_on="protein1", right_on="ensembl")
    pp    = pd.merge(pp, alias, left_on="protein2", right_on="ensembl")
    pp    = pp[["alias_x", "alias_y", "combined_score"]]
    
    pp.columns = ["gene_1", "gene_2", "score"]
    
    return(pp)

def string_target_binary_features(target_table, pp_table, th=900, output="table"):
    # Obtain binary matrix for gene target presence/absence

    pp_table     = pp_table.loc[pp_table.score>=th]
    # all_genes    = list(set(list(pp_table.gene_1) + list(pp_table.gene_2)))
    
    # target_table = target_table.loc[target_table.target.isin(all_genes)]
    
    all_drugs    = set(target_table.Compound)
    
    main_dict    = {}
    for i in all_drugs:

    	# Get target_genes
        target_genes  = list(set(target_table.loc[target_table.Compound==i].target))
        
        # Get partner genes
        partner_genes = pp_table.loc[(pp_table.gene_1.isin(target_genes)) | (pp_table.gene_2.isin(target_genes))]
        partner_genes = list(set(list(partner_genes.gene_1) + list(partner_genes.gene_2)))
        
        main_dict.update({i: list(set(target_genes+partner_genes))})
    
    new_set          = pd.DataFrame([[j,i] for j in main_dict.keys() for i in main_dict[j]], columns=["Compound", "target"])
    new_set["value"] = 1
    
    if output=="table":
        return(new_set) #[Compound, target, value==1]
    elif output=="pivot":
        return(new_set.pivot_table(index="Compound", columns="target", values="value", fill_value=0))

def target_lobico_parse(lobico_table, target_table):
	# Preps target table for compatibility with lobico-like table [cell_name, Compound, value]

	drugs   	 = list(set(lobico_table.Compound).intersection(target_table.index))
	lobico_table = lobico_table.loc[lobico_table.Compound.isin(drugs)]
	target_table = target_table.loc[drugs]

	return lobico_table, target_table

def string_target_expression_features(binary_target, exp_file, lobico_proc, output="table"):
    # exp_file such as 070818_cgp_exp.txt
    # binary_target: takes output (table) from string_target_binary_features()
    # lobico_proc: takes output from drug_exp_lobico_parse_v2()
    # Produces: Compound target expression features
    
    exp_feat = pd.read_csv(exp_file, sep="\t") # samples x genes
    
    # Consolidate genes on expression and target table
    genes          = list(set(list(exp_feat.columns)).intersection(binary_target.target))
    exp_feat       = exp_feat[genes]
    binary_target  = binary_target.loc[binary_target.target.isin(genes)]
    
    # Consolidate cells on expression and activity (lobico) table
    cells          = set(lobico_proc.cell_name).intersection(exp_feat.index)
    exp_feat       = exp_feat.loc[cells]
    lobico_proc    = lobico_proc.loc[lobico_proc.cell_name.isin(cells)]
    lobico_proc.columns = ["cell_name", "Compound", "activity"] #Rename at the end
    
    # Transform expression matrix to table [cell_name, target, value]
    exp_feat["cell_name"] = exp_feat.index
    exp_feat       = pd.melt(exp_feat, id_vars="cell_name", var_name="target", value_name="value")
    
    # Merge compound target to cells [Compound, cell_name, target, activity, value]
    binary_target = pd.merge(binary_target[["Compound", "target"]], lobico_proc, on="Compound") #[Compound, target, cell_name, activity]
    binary_target = pd.merge(binary_target, exp_feat, on=["target", "cell_name"])
    
    # Pivot manually two dimentional table (per compound)
    drugs         = set(binary_target.Compound)
    main_list     = []
    for d in drugs:
        t = binary_target.loc[binary_target.Compound==d].pivot_table(values="value", 
                                                                 index="cell_name", columns="target", fill_value=0)
        main_list.append(t.reset_index().assign(Compound=d))
    
    main_list     = pd.concat(main_list).fillna(0) #[Compound, cell_name, genes(targetxvalue)]
    main_list     = pd.merge(main_list, lobico_proc, on=["cell_name", "Compound"])
    
    # Separate features and get lobico back
    lobico_proc   = main_list[["cell_name", "Compound", "activity"]]
    main_list     = main_list.drop(["cell_name", "Compound", "activity"], axis=1)
    
    lobico_proc.columns = ["cell_name", "Compound", "value"]
    return lobico_proc, main_list

def string_target_expression_features_all(binary_target, exp_file, lobico_proc, output="table"):
	# NOTE: Similar to string_target_expression_features() but all gene features across compounds are available
	# 		unlike string_target_expression_features(), where 0s are filled for non-target genes
    # exp_file such as 070818_cgp_exp.txt
    # binary_target: takes output (table) from string_target_binary_features()
    # lobico_proc: takes output from drug_exp_lobico_parse_v2()
    # Produces: Compound target expression features
    
    exp_feat = pd.read_csv(exp_file, sep="\t") # samples x genes
    
    # Consolidate genes on expression and target table, Filtering genes
    genes          = list(set(list(exp_feat.columns)).intersection(binary_target.target))
    exp_feat       = exp_feat[genes]
    binary_target  = binary_target.loc[binary_target.target.isin(genes)] #[Compound, target, value==1]
    
    # Consolidate cells on expression and activity (lobico) table
    cells          = set(lobico_proc.cell_name).intersection(exp_feat.index)
    exp_feat       = exp_feat.loc[cells]
    lobico_proc    = lobico_proc.loc[lobico_proc.cell_name.isin(cells)]
    
    # Return lobico and feature space
    return lobico_proc, exp_feat

def Function_drug_name_zhang_to_gdsc(drug_name, replacement="gtz"):

	# Clean up first
	clean_up         = {"Zibotentan, ZD4054":"Zibotentan", "PXD101, Belinostat": "Belinostat"}
	if drug_name in clean_up.keys():
		return(clean_up[drug_name])

	# Then build to rename
	gdsc_to_zhang    = {"AICA Ribonucleotide":"AICAR","BAY-61-3606":"BAY 61-3606", "BX795":"BX-795",
						"CCT-018159":"CCT018159", "EHT-1864": "EHT 1864", "GSK1904529A": "GSK-1904529A",
						"GW441756":"GW 441756", "HG6-64-1":"HG-6-64-1", "Nutlin-3a (-)":"Nutlin-3a",
						"PD0325901":"PD-0325901", "PD173074":"PD-173074", "PLX-4720":"PLX4720",
						"SB216763":"SB 216763", "SB505124":"SB-505124", "SL0101":"SL 0101-1",
						"XAV939":"XAV 939", "ZM447439":"ZM-447439",
						"Tanespimycin":"17-AAG", "Wee1 Inhibitor":"681640", "Navitoclax":"ABT-263",
						"Linifanib":"ABT-869", "Veliparib":"ABT-888", "Quizartinib":"AC220",
						"Rucaparib":"AG-014699", "Motesanib":"AMG-706", "Ponatinib":"AP-24534",
						"Tretinoin":"ATRA", "Luminespib":"AUY922", "Tivozanib":"AV-951", "Saracatinib":"AZD-0530",
						"Selumetinib":"AZD6244", "Talazoparib":"BMN-673", "Avagacestat":"BMS-708163",
						"Idelalisib":"CAL-101", "Lestaurtinib":"CEP-701", "Alectinib":"CH5424802",
						"Pelitinib":"EKB-569", "Selisistat":"EX-527",
						"Daporinad":"FK866", "Pictilisib":"GDC0941", "Omipalisib":"GSK2126458",
						"Serdemetan":"JNJ-26854165", "Dacinostat":"LAQ824", "Enzastaurin":"LY317615",
						"Pevonedistat":"MLN4924", "Amuvatinib":"MP470", "Entinostat":"MS-275",
						"Dactolisib":"NVP-BEZ235", "Linsitinib":"OSI-906", "Palbociclib":"PD-0332991",
						"Refametinib":"RDEA119", "Seliciclib":"Roscovitine", "Ispinesib Mesylate":"SB-715992",
						"Fedratinib":"TG101348", "Tozasertib":"VX-680", "Cabozantinib":"XL-184",
						"Foretinib":"XL-880", "Sepantronium bromide":"YM155"}

	zhang_to_gdsc    = {y:x for x,y in gdsc_to_zhang.iteritems()}

	if replacement=="gtz":
		# print("renaming GDSC compounds names to published Zhang 2018 names")
		if drug_name in gdsc_to_zhang.keys():
			# print("Found compatible drug name")
			return(gdsc_to_zhang[drug_name])
		else:
			# print("no proper name needed or found")
			return(drug_name)

	elif replacement=="ztg":
		# print("renaming Zhang published names to original GDSC names")
		if drug_name in zhang_to_gdsc.keys():
			# print("Found compatible drug name")
			return(zhang_to_gdsc[drug_name])
		else:
			# print("no proper name needed or found")
			return(drug_name)

def cgp_pubchem_to_smiles(file_in, file_in_updated):
	# file_in as in "drug_ids.csv" or "pubchem_id.csv" (actually same file, obtained form gdsc) 
	# NOTE: file_in can also be a cleaned version of pubchem csv table dowloaded straight from gdsc, we just need a pubchem ID
	# UPDATE: Added table downloaded straight from gdsc compound search (empty search): x_update
	# file_in_updated as in gdsc_drug_list.csv
    import pubchempy as pcp    
    
    # Load files
    x         = pd.read_csv(file_in, header=None, names=["Compound", "count", "ID","N", "DROP"]).drop_duplicates() # SOURCE 1
    x_update  = pd.read_csv(file_in_updated, usecols=[1,5], header=0, names=["Compound", "ID"]).drop_duplicates() # SOURCE 2

    # Clean up
    x["Compound"] = [i.rstrip("'") for i in list(x.Compound)]
    x.ID          = x.ID.astype("str")
    x_update      = x_update.query("ID!='none'").query("ID!='several'")
    x_update      = x_update.dropna()
    
    # Combine with updated
    x_update  = x_update.loc[~x_update.Compound.isin(list(x.Compound))]
    x         = pd.concat([x[["Compound", "ID"]], x_update])

    #NOTE: Manual clean up of IDs due to spotted technical duplicates#
    remove_ids = ["681640", "53298813", "57370134", "16760671", "447912", "16683866", "6851740"]
    x          = x.loc[~x.ID.isin(remove_ids)]
    x          = pd.concat(x,
    	pd.DataFrame([["T0901317", "447912"]], columns=x.columns)
    	)
    x          = x.drop_duplicates()
	##################################################################
    
    print(x)
    print(x.Compound.value_counts().reset_index())
    main_list = []
    for j in set(x.Compound):
        p_id = x.loc[x.Compound==j].ID.values[0]
        print(j, p_id)
        c = pcp.Compound.from_cid(p_id)
        main_list.append([ j.rstrip("'"),[int(i) for i in str(c.cactvs_fingerprint)]
                         ])
    
    main_list = pd.DataFrame([i[1] for i in main_list], 
                             index=[i[0] for i in main_list],
                            columns = ["s_%s"%j for j in xrange(len(main_list[0][1]))] )
    
    # Save as main_list.to_csv("...", sep=",", header=True, index_label="Compound")
    print("Done")
    return(main_list)

def keras_fc_autoencoder(input_feat, input_noisy, encode, decode, keepprob, slr, batch_size):
	# Fully connected autoencoder with early stopping (10%) and patience=10

	from keras.layers import Input, Dense, Dropout, BatchNormalization
	from keras.layers.advanced_activations import PReLU
	from sklearn.metrics import mean_squared_log_error,mean_squared_error, r2_score,mean_absolute_error
	from keras.optimizers import Nadam, Adam,SGD,Adagrad,Adadelta,RMSprop
	from keras.callbacks import ReduceLROnPlateau, LearningRateScheduler, EarlyStopping
	from keras.models import Model
	from keras import backend as K
	import gc

	print("Encoding layers")
	
	# Start
	input_shape = Input(shape=(input_feat.shape[1],))
	print(encode)
	e           = encode[0] # In case single layer
	vars()["encoded_%s"%e] = Dense(encode[0], activation=PReLU())(input_shape)
	vars()["e_batch_%s"%e] = BatchNormalization()(vars()["encoded_%s"%e])
	vars()["e_drop_%s"%e]  = Dropout(keepprob)(vars()["e_batch_%s"%e])

	# Decoding layers
	prev_layer = encode[0]
	print("prev_layer: ", prev_layer)
	for e in encode[1:]:
		print(e)
		vars()["encoded_%s"%e] = Dense(e, activation=PReLU())(vars()["e_drop_%s"%prev_layer])
		vars()["e_batch_%s"%e] = BatchNormalization()(vars()["encoded_%s"%e])
		vars()["e_drop_%s"%e]  = Dropout(keepprob)(vars()["e_batch_%s"%e])
		prev_layer = e

	# Encode layers
	vars()["d_drop_%s"%prev_layer] = vars()["e_drop_%s"%e]
	d = prev_layer  # In case single layer
	for d in decode[:-1]:
		print(d, vars()["d_drop_%s"%prev_layer])
		vars()["decoded_{}".format(d)]  = Dense(d, activation=PReLU())(vars()["d_drop_{}".format(prev_layer)])
		vars()["d_batch_%s"%d] = BatchNormalization()(vars()["decoded_%s"%d])
		vars()["d_drop_%s"%d]  = Dropout(keepprob)(vars()["d_batch_%s"%d])
		prev_layer = d

	vars()["last"] = Dense(decode[-1], activation="sigmoid")(vars()["d_drop_%s"%d])

	# Map the autoencoder
	autoencoder = Model(input_shape, vars()["last"])

	# Create encoding layer
	encoder     = Model(input_shape, vars()["encoded_%s"%e])

	# Optimizer
	optimizer = Nadam(lr=slr, beta_1=0.9, beta_2=0.999)

	# Compile
	autoencoder.compile(optimizer=optimizer, loss="mean_squared_error", metrics=["mse"])

	# Run model
	# reduce_lr      = LearningRateScheduler(lambda x: 1e-3 * 0.9 ** x)
	early_stopping = EarlyStopping(monitor='val_loss', patience=30)

	results        = autoencoder.fit(input_noisy, input_feat, 
								batch_size=batch_size, epochs=1000, validation_split=0.15, verbose=2,
								callbacks=[early_stopping])

	# Obtain prediction
	output          = encoder.predict(input_feat)
	K.clear_session()

	return output, autoencoder, results

def keras_fc_autoencoder_wrap(data_feat, arch, slr, batch_size, noise, keepprob):
	# Wrapper for keras_fc_autoencoder_wrap
	# Processes feature space and architecture

	import numpy as np

	init_features = data_feat.shape[1]
	print("init arch: ", init_features, arch)
	layers        = arch_layers(init_features, arch) #[i, i/2, i/4]
	rev_layers    = [i for i in reversed(layers)] #[i/4, i/2, i]

	train_noisy   = data_feat + noise * np.random.normal(loc=0.0, scale=1.0, size=data_feat.shape)
	train_noisy   = np.clip(train_noisy, 0., 1.) #Clip to range to 0-1

	autoencoded, model, results  = keras_fc_autoencoder(data_feat, train_noisy, layers[1:], rev_layers[1:], keepprob, slr, batch_size)

	return autoencoded, model, results

def parse_target_features_dae(lobico, target_features, exp_target, string_th):
    # Build feature space
    # Do we need to add target features??
    if target_features==True:
        print("Using drug target features")
        pp_hugo             = string_pp_to_hugo("/home/ubuntu/STRING/9606.protein.links.v10.5.txt",
                                                "STRING/9606.protein.aliases.v10.5.txt")
        drug_target         = process_drug_target("/home/ubuntu/CGP_FILES/Screened_Compounds.csv")

        if exp_target==True:
            print("Using target expression features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="table")
            lobico, target_feat = string_target_expression_features(drug_target, "CGP_FILES/070818_cgp_exp.txt", lobico)

            drugs_index         = list(target_feat.index)
            data_feat           = target_feat.reset_index(drop=True)
        else:
            print("Using target binary features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="pivot")
            lobico, target_feat = target_lobico_parse(lobico, drug_target)

            drugs_index         = list(target_feat.index)
            data_feat           = target_feat.reset_index(drop=True)
    
    lobico        = lobico.reset_index(drop=True)
    print(data_feat.head())
    print("All data:", data_feat.shape)
    return lobico, data_feat, drugs_index

def load_data(exp_arch, exp_gfilter, drug_arch, problem="regression", in_folder="/home/ubuntu/"):
    # problem: regression, classification
    # Data for ROC
    lobico_ori = process_paper_lobico("{}PAPERS/LORIO/lobico.csv".format(in_folder))
    lobico_ori["Compound"] = [Function_drug_name_zhang_to_gdsc(i, "ztg") for i in lobico_ori["Compound"]]

    ####### Prep data #####
    # Expression features
    exp_data      = exp_encode_parse("c2setcover", exp_gfilter, exp_arch, 0.5, in_folder=in_folder)
    
    # Mutation features
    mut_data      = process_gdsc_variant("{}CGP_FILES/gdsc_WES_variants.csv".format(in_folder), output="matrix", th=5)
    
    # Drug features
    if drug_arch=="original":
        drug_data     = pd.read_csv("{}PAPERS/ZHANG_2018/drug_feats_proc.txt".format(in_folder), sep="\t", index_col=0)
        
    elif drug_arch=="pubchem_smiles":
        drug_data     = pd.read_csv("{}CGP_FILES/pubchem_smiles.txt".format(in_folder), index_col=0)
        
    elif drug_arch=="original_pubchem_smiles":
        print("using original and pubchem smiles")
        original  = pd.read_csv("{}PAPERS/ZHANG_2018/drug_feats_proc.txt".format(in_folder), sep="\t", index_col=0)
        pubchem   = pd.read_csv("{}CGP_FILES/pubchem_smiles.txt".format(in_folder), index_col=0)
        drug_data = pd.merge(original, pubchem, left_index=True, right_index=True)
        
    elif "_arch_" in drug_arch:
        print("using combined datasets")
        c_list    = [i for i in drug_arch.split("_arch_")]
        pubchem   = pd.read_csv("{}CGP_FILES/pubchem_smiles.txt".format(in_folder), index_col=0)
        encoded   = drug_encode_parse(c_list[1], 0.5, in_folder=in_folder)
        drug_data = pd.merge(pubchem, encoded, left_index=True, right_index=True)
    else:
        drug_data     = drug_encode_parse(drug_arch, 0.5, in_folder=in_folder)

    # Is this for classification or regression? 
    if problem=="regression":
        lobico = cgp_act_post_process(pd.read_csv("{}CGP_FILES/v17.3_fitted_dose_response.csv".format(in_folder)),
                                        zscoring=False)
    elif problem=="classification":
        lobico = lobico_ori

    lobico, exp_data, mut_data, drug_data = drug_exp_lobico_parse_v2(lobico, exp_data, mut_data, drug_data, "ztg")

    return lobico, exp_data, mut_data, drug_data, lobico_ori

def load_data_v2(cell_sources, drug_sources, problem="regression", in_folder="/home/ubuntu/"):
	# problem: regression, classification
	# Updated to load all data needed for problem

	# Load cell expression sources first
	print("Cell sources used :")
	cell_features = []
	for i in cell_sources:
		if i=="raw":
			print(i)
			cell_features.append(
				pd.read_csv("{}CGP_FILES/070818_cgp_exp.txt".format(in_folder), sep="\t")
				)

		elif i=="gsea_dae":
			print(i)
			exp_gfilter = int(10) # PRE-DEFINED
			exp_arch    = "2"	  # PRE-DEFINED
			print("Using as gsea filter {} and DAE archictecture {}".format(exp_gfilter, exp_arch))
			cell_features.append(
				exp_encode_parse("c2setcover", exp_gfilter, exp_arch, 0.5, in_folder=in_folder)
				)

		elif i=="gdsc_var":
			print(i)
			cell_features.append(
				process_gdsc_variant("{}CGP_FILES/gdsc_WES_variants.csv".format(in_folder), output="matrix", th=5)
				)

		elif i=="cgc_giant":
			print("Expression features filtered by Cancer Census Genes thresholded at GIANT interactions")
			th = 0.98
			print("GIANT thresholded at {}".format(th))
			genes = process_cgc_file("{}CGC_FILES/cancer_gene_census.csv".format(in_folder))[1]
			cell_features.append(
				cgc_gene_giant_features(genes, th=th)
				)
		else:
			print("cell feature type '{}' not found".format(i))

	# Load drug features next
	print("Drug sources used :")
	drug_features = []
	for j in drug_sources:
		if j=="zhang_original":
			print("Features as used in Zhang et al. (PADEL features)")
			drug_features.append(
				pd.read_csv("{}PAPERS/ZHANG_2018/drug_feats_proc.txt".format(in_folder), sep="\t", index_col=0)
				)

		elif j=="zhang_dae":
			drug_arch = "8_16"
			print("Using DAE encoded features with architecture {}".format(drug_arch))
			drug_features.append(
				drug_encode_parse(drug_arch, 0.5, in_folder=in_folder)
				)

		elif j=="pubchem_smiles":
			print(j, "update pubchem")
			# pubchem_feat = pd.read_csv("{}CGP_FILES/pubchem_smiles.txt".format(in_folder), index_col=0)
			# pubchem_feat = pd.read_csv("{}CGP_FILES/pubchem_smiles_updated.txt".format(in_folder), index_col=0)
			# pubchem_feat.index = [Function_drug_name_zhang_to_gdsc(i, "gtz") for i in list(pubchem_feat.index)]
			# pubchem_feat = pubchem_feat.rename_axis("Compound")
			# drug_features.append(pubchem_feat)

			train_smiles   = load_gdsc_smiles(in_folder)
			train_fp       = construct_fingerprint_features(train_smiles)
			drug_features.append(train_fp)

		elif j=="target_string":
			print("Using GDSC drug-target combined with STRING connectivity dataset as binary features")
			string_th   = 900
			print("Using STRING cutoff at {}".format(string_th/1000.0))

			pp_hugo     = string_pp_to_hugo("{}STRING/9606.protein.links.v10.5.txt".format(in_folder),
											"{}STRING/9606.protein.aliases.v10.5.txt".format(in_folder))
			drug_target = process_drug_target("{}CGP_FILES/Screened_Compounds.csv".format(in_folder))
			drug_target = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="pivot")
			drug_target = drug_gdsc_zhang_rename_index(drug_target, "gtz")

			drug_features.append(drug_target)

		elif j=="target_string_dae":
			print("Using DAE encoded GDSC drug-target combined with STRING connectivity dataset")
			string_th   = 550
			encoded_target_arch = "2_16"
			print("Using STRING cutoff at {} and DAE architecture of {}".format( (string_th/1000.0), encoded_target_arch))

			target_feat = pd.read_csv("{}CGP_FILES/drug_target_encoded_binary_{}_{}.txt".
									  format(in_folder, string_th, encoded_target_arch), 
									  sep="\t", index_col=0)
			target_feat = drug_gdsc_zhang_rename_index(target_feat, "gtz")
			drug_features.append(target_feat)

		elif j=="target_giant_dae":
			print("Using DAE encoded GDSC drug-target combined with GIANT connectivity dataset")
			giant_th    = 0.7
			encoded_target_arch = "4_16"
			print("Using GIANT cutoff at {} and DAE architecture of {}".format(giant_th, encoded_target_arch))

			giant       = pd.read_csv("{}GIANT_FILES/giant_drug_target_encoded_{}_{}.txt".
									  format(in_folder, giant_th, encoded_target_arch), 
									  sep="\t", index_col=0)
			giant       = drug_gdsc_zhang_rename_index(giant, "gtz")
			drug_features.append(giant)

		else:
			print("drug feature type '{}' not found".format(j))

	# Obtain combined feature sources
	cell_features = pd.concat(cell_features, axis=1)
	cell_features = cell_features.dropna()

	for i in drug_features:
		print(i.shape)
		print(i.drop_duplicates().shape)
		print(len(set(i.index)))
	drug_features = pd.concat(drug_features, axis=1)
	drug_features = drug_features.dropna()

	# Finally lobico and target files
	lobico_ori = process_paper_lobico("{}PAPERS/LORIO/lobico.csv".format(in_folder))
	lobico_ori["Compound"] = [Function_drug_name_zhang_to_gdsc(i, "ztg") for i in lobico_ori["Compound"]]

	# Is this for classification or regression?
	if problem=="regression":
		lobico = cgp_act_post_process(pd.read_csv("{}CGP_FILES/v17.3_fitted_dose_response.csv".format(in_folder)),
									  zscoring=True)
	elif problem=="classification":
		lobico = lobico_ori

	# Lastly, clean up correct indexes and names for features and target variables
	print("Number of cell expression features: {}".format(cell_features.shape[1]))
	print("Number of drug molecular features: {}".format(drug_features.shape[1]))
	lobico, cell_features, drug_features = drug_exp_lobico_parse_v3(lobico, cell_features, drug_features, "ztg")

	return lobico, cell_features, drug_features, lobico_ori

def parse_features_v2(lobico, cell_feat, drug_feat=None):
	# Takes input from load_data_v2()

	if drug_feat is not None:
		data_feat = pd.concat([cell_feat.loc[lobico.cell_name,].reset_index(drop=True),
							   drug_feat.loc[lobico.Compound,].reset_index(drop=True)], 
							   axis=1)
	else:
		print("No drug features used!!!")
		data_feat = cell_feat.loc[lobico.cell_name,].reset_index(drop=True)

	lobico    = lobico.reset_index(drop=True)

	# print(data_feat.head())
	print("All data:", data_feat.shape)

	return lobico, data_feat

def process_gdsc_variant(gdsc_file, output="table", th=10):
    # gdsc file such as gdsc_WES_variants.csv

    gdsc = pd.read_csv(gdsc_file, usecols=[0,2,3,5,6,7])
    gdsc["effect_name"] = gdsc["Gene"] + "_" + gdsc["AA"]
    gdsc = (gdsc.effect_name.value_counts().reset_index().
            query('effect_name>=%s'%th).
            pipe(lambda x: gdsc.loc[gdsc.effect_name.isin(x["index"])])
           )
    
    if output=="table":
        return(gdsc)
    elif output=="matrix":
        gdsc["value"] = 1
        gdsc = gdsc.pivot_table(index="SAMPLE", columns="effect_name", values="value", fill_value=0)
        return(gdsc)

def keras_auc(y_true, y_pred):
    auc = tf.metrics.auc(y_true, y_pred)[1]
    K.get_session().run(tf.local_variables_initializer())
    return auc

def prep_output_files_keras(out_folder, out_var, exp_arch, exp_gfilter, drug_arch, keepprob, arch, slr, 
                                 target_features, exp_target, encoded_target, encoded_target_arch, string_th, gdsc_variants, giant_features,
                                 problem, drug):
	# out_var: label of file purpose (ie: drugscv, drugout)
    # Prep output files
    pred_file     = "{}cgpmlp_exp_drug_reg_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_txt".format(out_folder, out_var, exp_arch, exp_gfilter, 
                                               drug_arch, keepprob, arch, slr, target_features, exp_target, 
                                               encoded_target, encoded_target_arch, string_th,gdsc_variants, giant_features, problem, drug)
    if problem=="classification":
        with open(pred_file, "w") as f:
            f.write("Compound\tcell_name\ttrue_value\tprob_0\tprob_1\tbest_metric")
    elif problem=="regression":
        with open(pred_file, "w") as f:
            f.write("Compound\tcell_name\ttrue_value\tprediction\tbest_metric")

    return pred_file

def prep_output_files_keras_v2(out_folder, out_var, cell_feat, drug_feat, keepprob, arch, slr, problem, drug, storage="txt"):
	# out_var: label of file purpose (ie: drugscv, drugout)
	# cell_feat and drug_feat are list of features
	# storage: txt, pickle
	cell_feat     = ("_").join(cell_feat)
	drug_feat     = ("_").join(drug_feat)

	# Prep output files
	pred_file     = "{}cgpmlp_exp_drug_{}_CELL_{}_DRUG_{}_PARAM_{}_{}_{}_{}_{}.{}".format(
		out_folder, out_var, cell_feat, drug_feat, keepprob, arch, slr, problem, drug, storage)

	if storage=="txt":
		if problem=="classification":
			with open(pred_file, "w") as f:
				f.write("Compound\tcell_name\ttrue_value\tprob_0\tprob_1\tbest_train\tbest_metric")
		elif problem=="regression":
			with open(pred_file, "w") as f:
				f.write("Compound\tcell_name\ttrue_value\tprediction\tbest_train\tbest_metric")

	return pred_file

def fix_cell_drug_feat_indices(feat_main, exp_feat, drug_feat, target_feat):
	# Constructs feature space for exp, drug and target features; and consolidates names in target table

    common_cells = set(list(feat_main.cell_name)).intersection(exp_feat.index)
    common_drugs = set(list(feat_main.Compound)).intersection(drug_feat.index)
    common_drugs = set(list(feat_main.Compound)).intersection(target_feat.index)
    feat_main    = feat_main.loc[feat_main.Compound.isin(common_drugs)]

    data_feat    = pd.concat([exp_feat.loc[feat_main.cell_name,].reset_index(drop=True),
                              drug_feat.loc[feat_main.Compound,].reset_index(drop=True),
                              target_feat.loc[feat_main.Compound,].reset_index(drop=True)], axis=1)
    return feat_main, data_feat

def parse_features(lobico, exp_data, mut_data, drug_data, target_features, exp_target, 
					encoded_target, encoded_target_arch, string_th, gdsc_variants, giant_features,
					in_folder="/home/ubuntu/"):
    # Build feature space
    
    # Build cell line feature space first
    if gdsc_variants==True:
        print("Using gdsc variant features")
        cell_data = pd.concat([exp_data, mut_data.loc[exp_data.index,]], 
                              axis=1)
    else:
        print("Not using gdsc variant features")
        cell_data = exp_data
    print("cell features: ", cell_data.shape[1])
        
    # Do we need to add target features??
    if target_features==True:
        print("Using drug target features")
        pp_hugo             = string_pp_to_hugo("{}STRING/9606.protein.links.v10.5.txt".format(in_folder),
                                                "{}STRING/9606.protein.aliases.v10.5.txt".format(in_folder))
        drug_target         = process_drug_target("{}CGP_FILES/Screened_Compounds.csv".format(in_folder))

        if exp_target==True:
            print("Using target expression features")
            drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="table")
            lobico, target_feat = string_target_expression_features(drug_target, "{}CGP_FILES/070818_cgp_exp.txt".format(in_folder), lobico)

            data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                        drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                        target_feat.reset_index(drop=True)], axis=1)
        else:
            print("Using target binary features")
            
            if encoded_target==True:
                print("Using encoded features")
                target_feat       = pd.read_csv("{}CGP_FILES/drug_target_encoded_binary_{}_{}.txt".format(in_folder, string_th, encoded_target_arch), 
                                                sep="\t", index_col=0)
                lobico, data_feat = fix_cell_drug_feat_indices(lobico, cell_data, drug_data, target_feat)

            else:
                print("Using non-encoded features")
                drug_target         = string_target_binary_features(drug_target, pp_hugo, th=string_th, output="pivot")
                lobico, target_feat = target_lobico_parse(lobico, drug_target)

                data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                            drug_data.loc[lobico.Compound,].reset_index(drop=True),
                                            target_feat.loc[lobico.Compound].reset_index(drop=True)], axis=1)
        
    else:
        print("No target features used")
        data_feat      = pd.concat([cell_data.loc[lobico.cell_name,].reset_index(drop=True),
                                    drug_data.loc[lobico.Compound,].reset_index(drop=True)], axis=1)
    
    # Do we need to add GIANT features
    if giant_features==True:
    	print("Using GIANT encoded features")
    	giant         = pd.read_csv("{}GIANT_FILES/giant_drug_target_encoded_{}_{}.txt".format(in_folder,0.7,"4_16"), sep="\t", index_col=0)
    	data_feat     = pd.concat([data_feat, 
    		giant.loc[lobico.Compound,].reset_index(drop=True)], axis=1)

    lobico        = lobico.reset_index(drop=True)
    print(data_feat.head())
    print("All data:", data_feat.shape)

    return lobico, data_feat

def process_hgnc_entrez_file(file_in):
    # File such as
    x = pd.read_csv(file_in, sep="\t",usecols=["Approved Symbol","Status","Entrez Gene ID"],
                    dtype={"Approved Symbol":"str","Status":"str","Entrez Gene ID":"str"})
    x = x.loc[x.Status=="Approved"]
    x = x.dropna()
    
    x.columns = ["gene", "status","entrez"]
    return x[["entrez", "gene"]]

def process_giant_entrez(giant_file, hgnc_file):
    # Process giant file and extract hugo symbols
    # giant_file such as "all_tissues_top"
    # hgnc_file such as "genes_entrez.txt" (custom dowloaded from HGNC)
    
    # Obtain genes to process
    x = pd.read_csv(giant_file, sep="\t", header=None, names=["entrez_1", "entrez_2", "prob"],
                           dtype={"entrez_1":"str", "entrez_2":"str", "prob":"float"}, engine="c")
        
    print("processing gene names")
    y = process_hgnc_entrez_file(hgnc_file)

    # Final table
    x = pd.merge(x, y, left_on="entrez_1", right_on="entrez")
    x = pd.merge(x, y, left_on="entrez_2", right_on="entrez")
    
    # Clean up
    x = x[["gene_x", "gene_y", "prob"]]
    x.columns = ["gene_1", "gene_2", "prob"]
    return(x)

def giant_th_analysis(giant):
    
    bins = [round(i,2) for i in np.arange(0.1, 1.0, 0.1)] + [round(i,2) for i in np.arange(0.95, 0.99, 0.01)]
    
    th_table = []
    for i in bins:
        i_len =(giant.loc[giant.prob>=i]
                .pipe(lambda p: len(set(p["gene_1"]) | set(p["gene_2"]))))
#         print(i, i_len)
        th_table.append([i, i_len])
    
    th_table = pd.DataFrame(th_table, columns=["th", "n_genes"])
    return th_table

def giant_gene_features(giant, drug_target, th=0.9):
    # drug_target: [Compound, target]
    # giant: [gene_1, gene_2, prob]
    # Build feature matrix using giant gene-gene interaction network based on threshold th
    # NOTE: Updated to integrate separatedly as two sources of features "Target" and "Partner" genes

    # Apply threshold
    giant        = giant.loc[giant.prob>=th]
    
    all_drugs    = set(drug_target.Compound)
    
    # First get target features
    drug_target["value"] = 1
    target_table = drug_target.pivot_table(index="Compound", columns="target", values="value", fill_value=0)
    target_table.columns = [i+"_T" for i in list(target_table.columns)]
    
    # Second, process partner features
    partner_table = {}
    for i in all_drugs:
    
        # Get target_genes
        target_genes  = list(set(drug_target.loc[drug_target.Compound==i].target))
        
        # Get partner genes
        partner_genes = giant.loc[(giant.gene_1.isin(target_genes)) | (giant.gene_2.isin(target_genes))]
        partner_genes = list(set(list(partner_genes.gene_1) + list(partner_genes.gene_2)))
        
        # Update
        partner_table.update({i: partner_genes})
    
    partner_table          = pd.DataFrame([[j,i] for j in partner_table.keys() for i in partner_table[j]], columns=["Compound", "target"])
    partner_table["value"] = 1
    partner_table          = partner_table.pivot_table(index="Compound", columns="target", values="value", fill_value=0)
    partner_table.columns  = [i+"_P" for i in list(partner_table.columns)]
    
    # Join tables
    main_table = pd.concat([target_table.loc[all_drugs],
                            partner_table.loc[all_drugs]], axis=1)
    
    # Add zeroes where necessary (
    # NOTE: Probably for compounds that have target information, but no GIANT target-partners
    main_table = main_table.fillna(0)

    print("Giant features: ", main_table.shape)
    return(main_table)

def encoded_giant_features(th=0.9, arch="9_16", in_folder="/home/ubuntu/"):
    import gc
    
    print("loading files...")
    drug_target     = process_drug_target("{}CGP_FILES/Screened_Compounds.csv".format(in_folder))
    giant           = process_giant_entrez("{}GIANT_FILES/all_tissues_top".format(in_folder), 
                                           "{}HGNC_FILES/genes_entrez.txt".format(in_folder))
    
    giant_table     = giant_gene_features(giant, drug_target, th)
    
    print("modeling...")
    slr             = 0.0001
    batch_size      = 20
    noise           = 0.5
    keepprob        = float(0.5)
    autoencoded, model, results  = keras_fc_autoencoder_wrap(giant_table, arch, slr, 
                                                             batch_size, noise, keepprob)
    
    print("encoding")
    encoded         = pd.DataFrame(autoencoded, index=giant_table.index)
    
    del autoencoded
    del model
    del results
    gc.collect()
    
    return encoded

def process_cgc_file(file_in):
    # file_in as in "cancer_gene_census.csv" from the CGC
    # output: table and unique gene list (includes missense mutations, deletions, translations)
    
    x = pd.read_csv(file_in, sep=",")
    y = list(set(x["Gene Symbol"].values))
    return x, y

def cgc_gene_giant_features(genes, th=0.99, in_folder="/home/ubuntu/"):
    # Build reduced gene expression feature space using giant-network connectivity
    # giant: [gene_1, gene_2, prob]
    
    # Filter by threshold
    giant    = process_giant_entrez("{}GIANT_FILES/all_tissues_top".format(in_folder), 
                                  "{}HGNC_FILES/genes_entrez.txt".format(in_folder))
    giant    = giant.loc[giant.prob>=th]
    
    # Filter by cgc genes
    giant    = giant.loc[(giant.gene_1.isin(genes)) | (giant.gene_2.isin(genes))]
    genes    = list(set( list(giant.gene_1) + list(giant.gene_2)))
    
    # Apply expression values
    gene_exp = pd.read_csv("{}CGP_FILES/070818_cgp_exp.txt".format(in_folder), sep="\t")
    genes    = list(set(genes).intersection(list(gene_exp.columns.values)))
    gene_exp = gene_exp[genes]
    
    # Return
    return gene_exp

def process_drugbank_xml(file_in):
	# Function to process xml file from DrugBank and obtain the names of compounds (and their synonyms)
	#  and their gene targets
    #file_in such as "DRUGBANK_FILES/full_database.xml"
    import xml.etree.ElementTree as et
    import pandas as pd
    
    # Load files
    db_xml = et.parse(file_in)
    root   = db_xml.getroot()
    
    # Process
    namespace = '{http://www.drugbank.ca}'
    main_list = []
    for i in root.getchildren():

        # Get name
        name   = [i.find("{}name".format(namespace)).text.rstrip(" ").lstrip(" ")]

        # Add synonyms
        syns   = i.find("{}synonyms".format(namespace))
        for s in syns.getchildren():
            name.append(s.text)

        # Get targets
        genes  = []
        for j in i.find("{}targets".format(namespace)):
            targets = j.findall(".//{}gene-name".format(namespace))
            genes   = genes + targets

        # Append if found
        if len(genes)>0:
            genes = [g.text for g in genes]

            main_list.append(pd.DataFrame([[n,g] for n in name for g in genes], columns=["Compound", "target"]))

    main_list = pd.concat(main_list)
    
    # Return
    return main_list

def process_stitch_target(drug_filter, stitch_target, stitch_alias, string_gene_alias):
	# drug_filter: list of drugs to filter for
	# stitch_aliases: aliases file from STITCH such as STITCH_FILES/chemical.aliases.v5.0.tsv (Takes long to load) 
	# stitch_target: drug-arget file from STICH  such as STITCH_FILES/9606.protein_chemical.links.v5.0.tsv
	# string_gene_alias: STRING alias file susch as STRING_FILES/9606.protein.aliases.v10.5.txt

	print("Loading files...")
	stitch_alias  = pd.read_csv(stitch_alias, sep="\t")
	stitch_target = pd.read_csv(stitch_target, sep="\t")

	print("Processing")
	common        = set(drug_filter).intersection(stitch_alias.alias.str.upper())

	# Filter stitch alias for drugs of interest
	stitch_alias  = stitch_alias.loc[stitch_alias.alias.str.upper().isin(common)]

	# Filter target for drugs of interest using alias names
	stitch_target = stitch_target.loc[stitch_target.chemical.isin(set(stitch_alias.flat_chemical).union(stitch_alias.stereo_chemical))]  

	# Get hugo identifiers for ensembl ids
	stitch_hugo   = string_to_hugo(list(set(stitch_target.protein)), string_gene_alias)

	# Apply hugo identifiers to stitch
	stitch_target = pd.merge(stitch_target, stitch_hugo, left_on="protein", right_on="ensembl")

	# Combine with drug alias names
	stitch_target = pd.merge(stitch_target,
							 pd.melt(stitch_alias, id_vars="alias", value_vars=["flat_chemical", "stereo_chemical"]),
							 left_on="chemical", right_on="value")

	# Clean up
	stitch_target = stitch_target[["alias_y", "alias_x", "combined_score"]]
	stitch_target = stitch_target.drop_duplicates()
	stitch_target.columns = ["Compound", "target", "score"]

	# Return
	return stitch_target

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

def construct_fingerprint_features(test_drug):
	# Constructs feature matrix out of 2-column pandas table of names and smiles

	all_smiles = list(test_drug.smiles)
	all_drugs  = list(test_drug.Compound)
	features   = change_smiles_morgan_fingerprint(all_smiles)

	feat_table = pd.DataFrame(features, index=all_drugs,
							  columns = ["s_%s"%j for j in xrange(len(features[0]))] )

	return feat_table

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

# CTRP-related Functions
def process_ctrp_drug_smiles_target(file_in):
    # file_in such as "v20.meta.per_compound.txt"
    x = pd.read_csv(file_in, sep="\t")
    x = x[["cpd_name", "gene_symbol_of_protein_target", "cpd_smiles"]]
    
    # Clean up
    x.columns = ["Compound", "target", "smiles"]
    x         = x.fillna("None")
    
    # Obtain targets
    y = pd.concat([pd.DataFrame({"Compound":i,
                                 "target":j.target.str.split(";").values[0]}) for i,j in x.groupby("Compound")])
    
    
    return x[["Compound", "smiles"]], y[["Compound", "target"]]

    # Written to:
    # - .../CTRP_FILES/drug_smiles.txt
    # - .../CTRP_FILES/drug_target.txt

# Cheat-Sheet
# df[df['model'].str.contains('ac')]
# .query("index.str.contains('s_')") #'query' in this case being a column name, and this step being in the process of a pipe
#  - NOTE, you might need .query("index.str.contains('s_')", engine="python") in some instances