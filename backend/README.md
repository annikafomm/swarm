- networks/     contains one GENIE3 network and one SPONGE network 
- datasets/     contains the xenium dataset as zarr file and a visium dataset

(networks/ and datasets/ might be too large)

- gene_set_ES/  contains a R project with two main files
                - genie_regulonES.Rmd --> calculates the VIPER score (and trying other scores)
                - spongEffects.Rmd --> experimenting with spongEffects
- calc_scores/  contains python scripts and outputs to calculate squidpy scores
                - calc_scores.py --> main script to calculate the squidpy scores
                - other files (except for h5ad files) are not important
