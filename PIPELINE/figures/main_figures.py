import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cPickle as pickle
from sklearn import preprocessing
from sklearn.manifold import TSNE
from scipy.stats import spearmanr

# Minor functions
def gdsc_tumor_type(file_in, th):
    # file_in such as CGP_FILES/gdsc_WES_variants.csv 
    x  = (pd.read_csv(file_in, usecols=["SAMPLE", "Cancer Type"])
              .dropna()
              .drop_duplicates()
              .rename(columns={"Cancer Type":"cancer"}))

    th = x.cancer.value_counts().reset_index().query("cancer>=@th")
    th = list(th["index"])
    
    x  = x.loc[x.cancer.isin(th)]
    
    return (x)

############## COMPARE CORRELATION OF CELL LINE EXPRESSION TO ACTIVITY ##############
cgp_act = cgp_act_post_process(pd.read_csv("{}CGP_FILES/v17.3_fitted_dose_response.csv".format("/home/ubuntu/")),
                               zscoring=True)

# Expression correlation (Previously calculated based on Pairwise Spearman across all cells)
with open("CGP_FILES/cell_exp_cor_spearman.pickle", "rb") as handle:
    cell_gene_exp_cor = pickle.load(handle)

sns.set(font_scale=0.8)
sns.set_style('whitegrid')
sns.clustermap(cell_gene_exp_cor, figsize=[14,14])

# Activity correlation
sns.set(font_scale=0.8)
sns.set_style('whitegrid')
sns.clustermap(cgp_act.pivot_table(columns="cell_name", index="Compound", values="value").corr(method="pearson").fillna(0)
	, figsize=[14,14])

# Correlation of both to see if expression explains activity


####### COMPARE GDSC EXPRESSION VS GSEA-DAE COMPRESSION IN TERMS OF PHENOTYPE RETENTION ########
in_folder  ="/home/ubuntu/"
gsea_dae   = exp_encode_parse("c2setcover", 10, "2", 0.5, in_folder=in_folder)
cell_exp   = pd.read_csv("{}CGP_FILES/070818_cgp_exp.txt".format(in_folder), sep="\t")
tumor_info = gdsc_tumor_type("CGP_FILES/gdsc_WES_variants.csv", 40)

tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=500)
gsea_tsne     = tsne.fit_transform(gsea_dae.values)

tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=500)
exp_tsne      = tsne.fit_transform(cell_exp.values)

# Plots
sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.subplots(figsize=(12,12))
ax = sns.scatterplot(data=pd.merge(pd.DataFrame(exp_tsne, index=cell_exp.index, columns=["s1", "s2"]).reset_index(),
                          tumor_info, left_on="index", right_on="SAMPLE").query("cancer!='ALL'"),
           x="s1", y="s2", hue="cancer", s=200, alpha=0.4, legend=None)
ax.set(xlabel='First component', ylabel='Second component')

sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.subplots(figsize=(12,12))
ax=sns.scatterplot(data=pd.merge(pd.DataFrame(gsea_tsne, index=gsea_dae.index, columns=["s1", "s2"]).reset_index(),
                          tumor_info, left_on="index", right_on="SAMPLE").query("cancer!='ALL'"),
           x="s1", y="s2", hue="cancer", s=200, alpha=0.4)
ax.set(xlabel='First component', ylabel='Second component')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.).texts[0].set_text("Cancer type")

####### COMPARE GDSC EXPRESSION VS GSEA-DAE COMPRESSION IN TERMS OF CELL-TO-CELL CORRELATION ########
gsea_dae_cor = gsea_dae.transpose().corr(method="spearman")
cell_exp_cor = cell_exp.transpose().corr(method="spearman")

cell_feat_cor = pd.merge(pd.melt(gsea_dae_cor.reset_index(), id_vars="index"),
                         pd.melt(cell_exp_cor.reset_index(), id_vars="index"), 
                         on=["index", "variable"])

sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.subplots(figsize=(16,12))

ax = sns.scatterplot(data = cell_feat_cor.query("index!=variable"),
               x="value_x", y="value_y")
ax.set(xlabel='Cell-to-Cell Correlation of Autoencoded Features', 
       ylabel='Cell-to-Cell Correlation of Gene Expression Features')
print(spearmanr(cell_feat_cor.query("index!=variable").value_x, 
				cell_feat_cor.query("index!=variable").value_y))

######## COMPARE GSEA-DAE/CGC-GENES VS WHOLE EXPRESSION PERFORMANCE ON GEELEHER RIDGE REGRESSION #########
df = pd.read_csv("FIGURE_DATA/ridge_tests") # Obtained from ridge_tests.py
# Use 1_2 or 5_4
ge = "doc"
ge_dae = "doc_5_4_dae"

sns.set(font_scale=1.5)
sns.set_style('whitegrid')
plt.subplots(figsize=(12,12))

ax = sns.lineplot(data=pd.concat([df.query("source==@ge_dae"), df.query("source==@ge")]),
            x="fpr", y="tpr", hue="source", ci=None, style="kind")
legend = ax.legend()
legend.get_title()

# legend(loc='lower right')
plt.setp(ax.get_legend().get_texts(), fontsize='12')
plt.setp(ax.get_legend().get_title(), fontsize='22')



legend.texts[0].set_text("Feature type")
legend.texts[1].set_text("Modular")
legend.texts[2].set_text("Expression - {0:.2f}".format(df.query("source==@ge").roc.values[0]))
legend.texts[3].set_text("Modular source - AUROC")
legend.texts[4].set_text("DAE - {0:.2f}".format(df.query("(source==@ge_dae) & (kind=='DAE')").roc.values[0]))
legend.texts[5].set_text("DAE + CGC - {0:.2f}".format(df.query("(source==@ge_dae) & (kind=='ALL')").roc.values[0]))
legend.texts[6].set_text("CGC - {0:.2f}".format(df.query("(source==@ge_dae) & (kind=='CGC')").roc.values[0]))

plt.plot([0, 1], [0, 1], linewidth=2, color="black", linestyle="dashed")
