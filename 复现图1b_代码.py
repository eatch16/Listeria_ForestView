# Initial imports
import glob
import os
from collections import OrderedDict, defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shapefile
import skbio
import statsmodels.stats.multitest as smt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from scipy.stats import ttest_rel
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu  
from scipy.stats import shapiro
from statsmodels.stats.multitest import multipletests

# Folder paths
out_path = "/projects/leaph/yingxian/AMR_Listeria/Data/Diversity_Analysis/Out/"
blast_sum_path = "/projects/leaph/yingxian/AMR_Listeria/Data/Gene_Detection/Out/"
card_path = "/projects/leaph/yingxian/AMR_Listeria/AMR_Listeria_add/CARD/Data/Out/"
shared = "/projects/leaph/shared/project_data/listeria_US/"
data_path = "/projects/leaph/yingxian/AMR_Listeria/Data/Diversity_Analysis/In/"

# prepare input data
env = pd.read_csv(shared + "environ_all.csv", index_col=0)
genome = pd.read_csv(shared + "Listeria_genomes.csv", index_col=0)
pd.merge(genome, env, how="left", left_on="Sample ID", right_index=True).to_csv(
    data_path + "genomes_env.csv"
)

### function to get the proportion
def get_ARG_proportion(argMat, argSpec):
    argMat.set_index("id", inplace=True)
    argMat[argMat > 1] = 1
    argMat.reset_index(inplace=True)
    argSpec.rename(columns={"Isolate ID": "id"}, inplace=True)
    final_data = pd.merge(argMat, argSpec[["id", "Species"]], on=["id"], how="inner")
    print(argMat.shape, argSpec.shape, final_data.shape)
    req_cols = [_ for _ in final_data.columns if _ not in ["id", "Species"]]
    final_data_1 = final_data.groupby("Species", as_index=False)[req_cols].agg("mean")
    temp = final_data_1[req_cols].sum()
    keep_cols = temp[temp > 0].index
    final_data_1["Tot_sum"] = final_data_1[keep_cols].sum(axis=1)
    final_data_1.sort_values("Tot_sum", ascending=False, inplace=True)
    final_data_1.set_index("Species", inplace=True)
    final_data_2 = final_data_1[final_data_1.sum(axis=1) > 0]
    return final_data_2[keep_cols]

arg_mat_1 = pd.read_csv(blast_sum_path + "ARG_functional_matrix.csv")
arg_spec = pd.read_csv(shared + "/Listeria_genomes.csv")

data_1 = get_ARG_proportion(arg_mat_1, arg_spec)
data_1 = data_1.rename(columns={"lmo0919": "lin", "lmo1695": "mprF"})
## Heat map
sns.set(rc={"figure.figsize": (5, 5)})
my_palette = sns.light_palette("#D98880", as_cmap=True)
ax = sns.heatmap(
    data_1,
    cmap=my_palette,
    square=True,
    annot=False,
    fmt=".2f",
    annot_kws={"size": 5},
    cbar_kws={"label": "Fraction of ARGs among genomes"},
)
plt.title("Functional ARGs", fontsize=15)
ax.set(ylabel=None)
sns.diverging_palette(145, 300, s=60, as_cmap=True)
plt.savefig(
    out_path + "Heatmap_Fraction_Functional_ARG.pdf", bbox_inches="tight", dpi=600
)
