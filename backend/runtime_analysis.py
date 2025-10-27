import os
import asyncio
import scanpy as sc
from app import calculate_scores_helper
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt


nhood_perms = [100, 500, 1000, 5000, 10000, 50000, 100000]
autocorr_perms = [None, 1, 10, 50, 100, 500, 1000]

visium_files = "../backend/datasets_prepro_new"

"""
for file in os.listdir(visium_files):
    if file.startswith("GSM"):
        
        file_path = os.path.join(visium_files, file)
        print(file_path)
        adata = sc.read_h5ad(file_path)
        print(adata.n_vars)
        print(len(set(adata.obs["leiden"])))

        for i in range(len(nhood_perms)):
            
            if autocorr_perms[i] == None:
                python_params = ["-input", file_path, 
                            "-moranI",
                            "-gearyC",
                            "-nhood_enrichment", "-n_perms_nhood", f"{nhood_perms[i]}"
                            ]
            else:
                python_params = ["-input", file_path, 
                                "-moranI", "-n_perms_autocorr_mI", f"{autocorr_perms[i]}",
                                "-gearyC", "-n_perms_autocorr_gC", f"{autocorr_perms[i]}",
                                "-nhood_enrichment", "-n_perms_nhood", f"{nhood_perms[i]}"
                                ]
            R_params = []
            
            params_dict = {"params_python_script": python_params,
                        "params_R_script": R_params}
            
            asyncio.run(calculate_scores_helper(os.path.join("../backend/datasets_runtime", f'{file.replace(".h5ad", "")}_n_perms_{i}'), params_dict))
"""

score_files = Path("../backend/datasets_runtime")

moran_I_dict = {}
moran_padj_dict = {}
geary_C_dict = {}
geary_padj_dict = {}
nhood_count_dict = {}
nhood_zscore_dict = {}

for h5ad_file in score_files.rglob('*.h5ad'):
    parent_dir_name = h5ad_file.parents[1].name
    parts = parent_dir_name.split('_')
    m_code = parts[1]
    n_perms = int(parts[-1])

    label_autocorr = f"{m_code}_{autocorr_perms[n_perms]}"
    label_nhood = f"{m_code}_{nhood_perms[n_perms]}"

    print(m_code, n_perms, h5ad_file)

    adata = sc.read_h5ad(h5ad_file)

    moran_I_dict[label_autocorr] = adata.uns["moranI"]["I"]
    moran_padj_dict[label_autocorr] = adata.uns["moranI"]["pval_norm"]
    geary_C_dict[label_autocorr] = adata.uns["gearyC"]["C"]
    geary_padj_dict[label_autocorr] = adata.uns["gearyC"]["pval_norm"]
    nhood_count_dict[label_nhood] = adata.uns["leiden_nhood_enrichment"].get("count")
    nhood_zscore_dict[label_nhood] = adata.uns["leiden_nhood_enrichment"].get("zscore")

scores_dict = {"moran_I": moran_I_dict,
               "moran_padj": moran_padj_dict,
               "geary_C": geary_C_dict,
               "geary_padj": geary_padj_dict,
               "nhood_count": nhood_count_dict,
               "nhood_zscore": nhood_zscore_dict
               }

"""
for key in scores_dict.keys():
    score_df = pd.DataFrame(scores_dict.get(key))

    for m in ["M2", "M5", "M6", "M8", "M13"]:
        selected_cols = [f"{m}_{i}" for i in autocorr_perms]
        score_df_filt = score_df.loc[:, score_df.columns.str.startswith(m)][selected_cols].dropna()

        # Compute correlation matrix and Convert to square matrix
        corr_matrix = squareform(pdist(score_df_filt.T, metric='euclidean'))
        labels = score_df_filt.columns 
        corr_df = pd.DataFrame(corr_matrix, index=labels, columns=labels)

        plt.figure(figsize=(10, 8))
        sns.heatmap(corr_df, annot=True, cmap='viridis', square=True)
        plt.title(f'Pairwise Euclidean Distances Between {key} for dataset {m}')
        plt.tight_layout()
        plt.savefig(f"../backend/datasets_runtime/heatmaps/{key}_{m}_distances.png", dpi=300)
        plt.show()
"""


def pairwise_mae(df1, df2):
    # Assumes df1 and df2 have the same shape and aligned indices/columns
    return np.mean(np.abs(np.asarray(df1) - np.asarray(df2)))

for key in scores_dict.keys():
    print(key)
    for m in ["M2", "M5", "M6", "M8", "M13"]:
        filtered_dict = {k: v for k, v in scores_dict.get(key).items() if k.startswith(m)}
        
        if not key.startswith("nhood"):
            keys = [f"{m}_{i}" for i in autocorr_perms]
            n = len(keys)

            score_df = pd.DataFrame(filtered_dict)[keys].dropna()

            # Prepare empty matrix for MAE values
            mae_matrix = pd.DataFrame(np.zeros((n, n)), index=keys, columns=keys)
            # Compute pairwise MAE
            for i in range(n):
                for j in range(i, n):
                    mae_value = pairwise_mae(score_df[keys[i]], score_df[keys[j]])
                    mae_matrix.iloc[i, j] = mae_value
                    mae_matrix.iloc[j, i] = mae_value  # symmetric
        else:
            keys = [f"{m}_{i}" for i in nhood_perms]
            n = len(keys)

            mae_matrix = pd.DataFrame(np.zeros((n, n)), index=keys, columns=keys)
            # Compute pairwise MAE
            for i in range(n):
                for j in range(i, n):
                    mae_value = pairwise_mae(filtered_dict[keys[i]], filtered_dict[keys[j]])
                    mae_matrix.iloc[i, j] = mae_value
                    mae_matrix.iloc[j, i] = mae_value  # symmetric
           
        plt.figure(figsize=(10, 8))
        sns.heatmap(mae_matrix, annot=True, fmt=".2f", cmap='viridis', square=True, vmin=0, vmax=1)
        plt.title(f'Pairwise Mean Absolute Error of {key} across conditions in dataset {m}')
        plt.tight_layout()
        plt.savefig(f"../backend/datasets_runtime/heatmaps_MAE/{key}_{m}_distances.png", dpi=300)
        
        