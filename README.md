![MOPITAS logo](logo.png)

# MOPITAS MaPra

This is the Github-Repository of the MoPiTas Project in the Master Practical course of DaisyBioLab 2025.

## Dataset-Table

https://docs.google.com/spreadsheets/d/1RH_lB4OEUEgZMIYu59_RTwcPmkC-ciyApEnlhvPSJ7E/edit?usp=sharing 

## Frontend

Start the frontend development server:

```sh
cd frontend
npm start
```
<img width="1072" height="921" alt="Screenshot from 2025-08-14 14-33-33" src="https://github.com/user-attachments/assets/145f8ff9-9150-4cda-9d2e-c2eb8f3f7c35" />

# Format Requirements
The web-tool is working with the anndata-format. This is a tutorial on how counts and scores should be saved.

**Counts** - anndata.X
**Clusterings/Observations (per Spot)** - anndata.obs (all observation columns will be shown so make sure to clean your anndata before upload). Save the Leiden-clustering for leiden-based scores under “leiden”
**Scores (Geary C, Moran I, leiden-centrality, leiden-cooccurence)** - anndata.uns as DataFrames / Dicts
**DF:** “gearyC” cols: C, pval_norm
**DF:** “moranI” cols: I, pval_norm
**DF:** “leiden_centrality_scores” cols: degree_centrality, average_clustering, closeness_centrality
**DICT:** “leiden_co_occurrence”: dict with intervals: n_intervals - 1 ,occ: 3D array (n_clusters,n_clusters, n_intervals) (Make sure to select the right num_intervals in Upload form)
**DICT**: “leiden_nhood_enrichment” dict: counts: (n_clust, n_clust), zscore(n_clust, nclust)
**AUcell Scores**
**DF:** cols “Regulon”: Genesets , rows “Cell”: cell_names as barcodes
**DICT:** genesets/ regulons should be saved as a dict with regulon_name: tuple(gene1,...,gene_n)

## Tangram (if mapping precomputed)
## LIANA+ Scores (if precomputed)

