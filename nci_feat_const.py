import numpy as np
import pandas as pd
import gc

TABLES = "/tigress/zamalloa/TABLES/TCGA.TRAINING/"

x = pd.read_csv(TABLES +  "nci60_all_feat_table_scaled_round.2.csv", sep="\t", nrows=100)
cosmic = pd.read_csv("TABLES/TCGA.TRAINING/cosmic.genes.csv", sep="\t", header=None)

met = x.columns.values[:864]
genes = x.columns.values[864:]
genes = list(set(genes) & set(list(cosmic[0])))

filt_cols = list(met) + genes

del x
del cosmic
gc.collect()

df_test  = pd.read_csv(TABLES + "nci60_all_feat_table_scaled_round.2.csv", sep="\t", nrows=5)
float_cols = [c for c in df_test] # if df_test[c].dtype=="float64"]
float32_cols = {c:np.float32 for c in float_cols}
df  = pd.read_csv(TABLES  + "nci60_all_feat_table_scaled_round.2.csv", sep="\t", engine="c", dtype=float32_cols)
gc.collect()

df = df[filt_cols]
gc.collect()

with open(TABLES + "nci60_all_feat_table_scaled_round.3.csv", "w") as ff:
    df.iloc[:30000,:].to_csv(ff, header=True, sep="\t", index=False, float_format='%.3f')

for i in xrange(30000, 630000, 30000):
    print(i, i+30000)
    with open(TABLES + "nci60_all_feat_table_scaled_round.3.csv", "a") as ff:
        df.iloc[i:(i+30000),:].to_csv(ff, header=False, sep="\t", index=False, float_format='%.3f')

with open(TABLES + "nci60_all_feat_table_scaled_round.3.csv", "a") as ff:
    df.iloc[63000:,:].to_csv(ff, header=False, sep="\t", index=False, float_format='%.3f')

print("Done writing feature space")
