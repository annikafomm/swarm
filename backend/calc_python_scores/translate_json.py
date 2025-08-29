
def dict2params(param_dict):
    python_params = []
    network_params = []

    for key in param_dict.keys():
        match key:
            case "spatial_h5ad":
                python_params.append("-input")
                python_params.append(param_dict.get(key))
            case "spatial_normalization":
                if param_dict.get(key):
                    python_params.append("-normalize_st")
            case "spatial_filtering":
                if param_dict.get(key):
                    python_params.append("-filter_st")
            
            case "use_tangram": 
                if param_dict.get(key):
                    python_params.append("-tangram")
                    network_params.append("--tangram")
            case "single_cell_h5ad": 
                python_params.append("-sc_path")
                python_params.append(param_dict.get(key))
            case "singlecell_filtering": 
                if param_dict.get(key):
                    python_params.append("-filter_sc")
            case "singlecell_normalization": 
                if param_dict.get(key):
                    python_params.append("-normalize_sc")
            
            case "score_network":
                if param_dict.get(key):
                    python_params.append("-R_scores")
            
            case "score_liana_plus":
                if param_dict.get(key):
                    python_params.append("-liana")
            case "liana_composition_column": 
                python_params.append("-cell_comp_key")
                python_params.append(param_dict.get(key))
            case "liana_genie3_network":
                python_params.append("-grn") 
                python_params.append(param_dict.get(key))
            case "liana_pathway_network": 
                python_params.append("-pathway_net")
                python_params.append(param_dict.get(key))

            case "squidpy_moranI": 
                if param_dict.get(key):
                    python_params.append("-moranI")
            case "squidpy_moranI_n_perms": 
                python_params.append("-n_perms_autocorr_mI")
                python_params.append(param_dict.get(key))
            case "squidpy_moranI_two_tailed": 
                if param_dict.get(key):
                    python_params.append("-two_tailed_mI")
            case "squidpy_moranI_corr_method": 
                python_params.append("-corr_method_mI")
                python_params.append(param_dict.get(key))

            case "squidpy_gearyC": 
                if param_dict.get(key):
                    python_params.append("-gearyC")
            case "squidpy_gearyC_n_perms": 
                python_params.append("-n_perms_autocorr_gC")
                python_params.append(param_dict.get(key))
            case "squidpy_gearyC_two_tailed": 
                if param_dict.get(key):
                    python_params.append("-two_tailed_gC")
            case "squidpy_gearyC_corr_method": 
                python_params.append("-corr_method_gC")
                python_params.append(param_dict.get(key))

            case "squidpy_centrality_score": 
                if param_dict.get(key):
                    python_params.append("-centrality_scores")
            case "squidpy_centrality_score_cluster_key": 
                python_params.append("-cluster_cs")
                python_params.append(param_dict.get(key))

            case "squidpy_co_occurrence": 
                if param_dict.get(key):
                    python_params.append("-co_occurrence")
            case "squidpy_co_occurrence_cluster_key": 
                python_params.append("-cluster_co")
                python_params.append(param_dict.get(key))
            case "squidpy_co_occurrence_interval": 
                python_params.append("-interval")
                python_params.append(param_dict.get(key))
            case "squidpy_co_occurrence_n_splits": 
                python_params.append("-n_splits")
                python_params.append(param_dict.get(key))

            case "squidpy_neighborhood_enrichment": 
                if param_dict.get(key):
                    python_params.append("-nhood_enrichment")
            case "squidpy_neighborhood_enrichment_cluster_key": 
                python_params.append("-cluster_nhood")
                python_params.append(param_dict.get(key))
            case "squidpy_neighborhood_enrichment_library_key": 
                python_params.append("-library_key")
                python_params.append(param_dict.get(key))
            case "squidpy_neighborhood_enrichment_n_perms": 
                python_params.append("-n_perms_nhood")
                python_params.append(param_dict.get(key))

            case "alg_viper": 
                if param_dict.get(key):
                    network_params.append("--viper")
            case "alg_aucell": 
                if param_dict.get(key):
                    network_params.append("--aucell")
            case "alg_gsva": 
                if param_dict.get(key):
                    network_params.append("--gsva")
            case "alg_ssgsea": 
                if param_dict.get(key):
                    network_params.append("--ssgsea")

            case "sponge_networkanalysis": 
                network_params.append("--sponge_analysis")
                network_params.append(param_dict.get(key))
            case "sponge_networkinteractions": 
                network_params.append("--sponge_network")
                network_params.append(param_dict.get(key))
            case "net_m_score_threshold": 
                network_params.append("--mscor")
                network_params.append(param_dict.get(key))
            case "net_p_adjust": 
                network_params.append("--padj")
                network_params.append(param_dict.get(key))
            case "net_ensembl_id_col": 
                network_params.append("--ensembl_col")
                network_params.append(param_dict.get(key))
            case "net_feature_col": 
                network_params.append("--feature_col")
                network_params.append(param_dict.get(key))
            case "net_rna_types": 
                network_params.append("--RNA_types")
                network_params.append(param_dict.get(key))
            case "net_max_modules": 
                network_params.append("--max_modules")
                network_params.append(param_dict.get(key))

            case "genie3_network": 
                network_params.append("--genie_network")
                network_params.append(param_dict.get(key))
            case "genie3_top_n_weights": 
                network_params.append("--top_n")
                network_params.append(param_dict.get(key))
            case "genie3_n_regulatory_genes": 
                network_params.append("--k_reg_genes")
                network_params.append(param_dict.get(key))
            case "genie3_n_regulons": 
                network_params.append("--n_regulons")
                network_params.append(param_dict.get(key))

    return (python_params, network_params)