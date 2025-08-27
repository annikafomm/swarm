
def dict2params(param_dict):
    tangram_params = []
    liana_params = []
    python_params = []
    network_params = []

    for key in param_dict.keys():
        match key:
            case "normalization":
                python_params.append("-normalize")
            case "filteringSpatial":
                python_params.append("-filter_st")
            case "filteringSingleCell":
                python_params.append("-filter_st")
            case "files":
                files_dict = param_dict.get(key)
                for fkey in files_dict.keys():
                    if files_dict.get(fkey) != None:
                        match fkey:
                            case "spatialFile":
                                python_params.append("-input")
                                python_params.append(files_dict.get(fkey))
                                tangram_params.append("-sp_path")
                                tangram_params.append(files_dict.get(fkey))
                            case "singleCellFile":
                                tangram_params.append("-sc_path")
                                tangram_params.append(files_dict.get(fkey))
                            case "precomputedFile":
                                pass
                            case "spongeNetworkAnalysis":
                                network_params.append("--sponge_analysis")
                                network_params.append(files_dict.get(fkey))
                            case "spongeNetworkInteractions":
                                network_params.append("--sponge_network")
                                network_params.append(files_dict.get(fkey))
                            case "genieFile":
                                network_params.append("--genie_network")
                                network_params.append(files_dict.get(fkey))
        
                        
            case "tangram":
                pass
            case "scores":
                pass
            case "options":
                options_dict = param_dict.get(key)
                for okey in options_dict.keys():
                    match okey:
                        case "tangram":
                            pass
                        case "liana":
                            pass
                        case "sponge":
                            score_dict = options_dict.get(okey)
                            for skey in score_dict.keys():
                                match skey:
                                    case "mScoreThreshold":
                                        network_params.append("--mscor")
                                    case "pAdj":
                                        network_params.append("--padj")
                                    case "featureKey":
                                        network_params.append("--feature_col")
                                    case "ensemblIdKey":
                                        network_params.append("--ensembl_col")
                                    case "rnaTypes":
                                        network_params.append("--RNA_types")
                                    case "maxModules":
                                        network_params.append("--max_modules")

                                network_params.append(score_dict.get(skey))
                                
                                # scores missing

                        case "genie": # does this exist?
                            pass
                        case "squidpy":
                            score_dict = options_dict.get(okey)
                            for skey in score_dict.keys():
                                match skey:
                                    case "clusterKey": 
                                        python_params.append("-cluster")
                                    case "coOccurInterval": 
                                        python_params.append("-interval")
                                    case "coOccurNSplits": 
                                        python_params.append("-n_splits")
                                    case "neighLibraryKey": 
                                        python_params.append("-library_key")
                                    case "neighNPerms":
                                        python_params.append("-n_perms_nhood")
                                    case "nPermutations":
                                        python_params.append("-n_perms_autocorr")
                                    case "twoTailed": 
                                        python_params.append("-two_tailed")
                                    case "corrMethod":
                                        python_params.append("-corr_method")
                                
                                python_params.append(score_dict.get(skey))
                    
                                    
                                    
                            
                        
            
        
{
    "options": {
        "tangram": {
        "cellTypeKey": "cell_type"
        },
        "sponge": {
        "enrichmentMethod": "GSVA",
        },
        "squidpy": {
        "method": "MoransI",
        },
        "liana": {
        "ctCompositionKey": "tangram_ct_pred"
        }
    },

  "scores": [
    "LIANA+",
    "SPONGeffects",      # rename to net_scores / network_scores (combined with VIPER)
    "squidpy",
    "VIPER"
  ],

  "useDefaultLiana": False
}