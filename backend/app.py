#from flask import Flask, request, jsonify
import os
import subprocess
import random
import shutil
from datetime import datetime
import asyncio

def dict2params(param_dict):
    python_params = []
    network_params = []
    multiome_params = []
    multiome_params_py = []
    compute_R_scores = False

    for key in param_dict.keys():
        match key:
            case "spatial":
                for skey in param_dict.get(key).keys():
                    if param_dict.get(key).get(skey):
                        match skey:
                            case "normalization":
                                python_params.append("-normalize_st")
                            case "filtering":
                                python_params.append("-filter_st")
            case "files":
                for fkey in param_dict.get(key).keys():
                    if param_dict.get(key).get(fkey) is not None:
                        match fkey:
                            case "spatial_h5ad":
                                python_params.append("-input")
                                python_params.append(param_dict.get(key).get(fkey))
                            case "single_cell_h5ad":
                                python_params.append("-sc_path")
                                python_params.append(param_dict.get(key).get(fkey))
                            case "multiome_rds":
                                multiome_params.append("--multiome_rds")
                                multiome_params.append(param_dict.get(key).get(fkey))
                            case "fragments_tsv_gz":
                                multiome_params.append("--fragments_tsv_gz")
                                multiome_params.append(param_dict.get(key).get(fkey))
                            case "fragments_tsv_gz_tbi":
                                multiome_params.append("--fragments_tsv_gz_tbi")
                                multiome_params.append(param_dict.get(key).get(fkey))
                            case "genie3_network":
                                network_params.append("--genie_network")
                                network_params.append(param_dict.get(key).get(fkey))
                            case "sponge_networkanalysis":
                                network_params.append("--sponge_analysis")
                                network_params.append(param_dict.get(key).get(fkey))
                            case "sponge_networkinteractions":
                                network_params.append("--sponge_network")
                                network_params.append(param_dict.get(key).get(fkey))
                            case "liana_genie3_network":
                                python_params.append("-grn")
                                python_params.append(param_dict.get(key).get(fkey))
                            case "liana_pathway_network":
                                python_params.append("-pathway_net")
                                python_params.append(param_dict.get(key).get(fkey))
            case "tangram":
                for tkey in param_dict.get(key).keys():
                    if tkey == "use" and param_dict.get(key).get(tkey):
                        python_params.append("-tangram")
                        network_params.append("--tangram")
                    else:
                        if param_dict.get(key).get(tkey) is not None and param_dict.get(key).get(tkey):
                            match tkey:
                                case "filtering":
                                    python_params.append("-filter_sc")
                                case "normalization":
                                    python_params.append("-normalize_sc")
            case "multiome":
                for tkey in param_dict.get(key).keys():
                    if tkey == "use" and param_dict.get(key).get(tkey):
                        multiome_params.append("--multiome")
            case "genome":
                if param_dict.get(key) is not None:
                    multiome_params.append("--genome")
                    multiome_params.append(param_dict.get(key))
            case "scores":
                for skey in param_dict.get(key).keys():
                    if param_dict.get(key).get(skey):
                        match skey:
                            case "network":
                                python_params.append("-R_scores")
                                compute_R_scores = True
                            case "liana_plus":
                                python_params.append("-liana")
                            case "chromVar":
                                multiome_params.append("--chromvar")
                                multiome_params_py.append("-chromvar")
                                multiome_params_py.append("-moranI")
                                multiome_params_py.append("-gearyC")
                            # case "differential_motif_activity":
                            #     multiome_params.append("--differential_motif_activity")
                            #     multiome_params_py.append("-differential_motif_activity")
                            case "motif_enrichment":
                                multiome_params.append("--motif_enrichment")
                                multiome_params_py.append("-motif_enrichment")
                            case "footprinting":
                                multiome_params.append("--footprinting")
                                multiome_params_py.append("-footprinting")

            case "network":
                for nkey in param_dict.get(key).keys():
                    match nkey:
                        case "algorithms":
                            for akey in param_dict.get(key).get(nkey).keys():
                                if param_dict.get(key).get(nkey).get(akey):
                                    match akey:
                                        case "viper":
                                            network_params.append("--viper")
                                        case "aucell":
                                            network_params.append("--aucell")
                                        case "gsva":
                                            network_params.append("--gsva")
                                        case "ssgsea":
                                            network_params.append("--ssgsea")
                        case "sponge_params":
                            for skey in param_dict.get(key).get(nkey).keys():
                                if param_dict.get(key).get(nkey).get(skey) is not None:
                                    match skey:
                                        case "m_score_threshold":
                                            network_params.append("--mscor")
                                            network_params.append(str(float(param_dict.get(key).get(nkey).get(skey))))
                                        case "p_adjust":
                                            network_params.append("--padj")
                                            network_params.append(str(float(param_dict.get(key).get(nkey).get(skey))))
                                        case "ensembl_id_col":
                                            network_params.append("--ensembl_col")
                                            network_params.append(param_dict.get(key).get(nkey).get(skey))
                                            python_params.append("-ensembl_col")
                                            python_params.append(param_dict.get(key).get(nkey).get(skey))
                                        case "feature_col":
                                            network_params.append("--feature_col")
                                            network_params.append(param_dict.get(key).get(nkey).get(skey))
                                            python_params.append("-feature_col")
                                            python_params.append(param_dict.get(key).get(nkey).get(skey))
                                        case "rna_types":
                                            network_params.append("--RNA_types")
                                            network_params.append(param_dict.get(key).get(nkey).get(skey))
                                        case "max_modules":
                                            network_params.append("--max_modules")
                                            network_params.append(str(int(param_dict.get(key).get(nkey).get(skey))))
                        case "genie3_params":
                            for gkey in param_dict.get(key).get(nkey).keys():
                                if param_dict.get(key).get(nkey).get(gkey) is not None:
                                    match gkey:
                                        case "top_n_weights":
                                            network_params.append("--top_n")
                                            network_params.append(str(int(param_dict.get(key).get(nkey).get(gkey))))
                                        case "n_regulatory_genes":
                                            network_params.append("--k_reg_genes")
                                            network_params.append(str(int(param_dict.get(key).get(nkey).get(gkey))))
                                        case "n_regulons":
                                            network_params.append("--n_regulons")
                                            network_params.append(str(int(param_dict.get(key).get(nkey).get(gkey))))
            case "squidpy":
                for skey in param_dict.get(key).keys():
                    if (not skey.endswith("_params")) and (param_dict.get(key).get(skey)):
                        match skey:
                            case "moranI":
                                python_params.append("-moranI")
                                score_dict = param_dict.get(key).get(f"{skey}_params")
                                for param in score_dict.keys():
                                    if param == "n_perms" and score_dict.get(param) is not None:
                                        python_params.append("-n_perms_autocorr_mI")
                                        python_params.append(str(int(score_dict.get(param))))
                                    if param == "two_tailed" and score_dict.get(param):
                                        python_params.append("-two_tailed_mI")
                                    if param == "corr_method" and score_dict.get(param) is not None:
                                        python_params.append("-corr_method_mI")
                                        python_params.append(score_dict.get(param))
                            case "gearyC":
                                python_params.append("-gearyC")
                                score_dict = param_dict.get(key).get(f"{skey}_params")
                                for param in score_dict.keys():
                                    if param == "n_perms" and score_dict.get(param) is not None:
                                        python_params.append("-n_perms_autocorr_gC")
                                        python_params.append(str(int(score_dict.get(param))))
                                    if param == "two_tailed" and score_dict.get(param):
                                        python_params.append("-two_tailed_gC")
                                    if param == "corr_method" and score_dict.get(param) is not None:
                                        python_params.append("-corr_method_gC")
                                        python_params.append(score_dict.get(param))

                            case "centrality_score":
                                python_params.append("-centrality_scores")
                                score_dict = param_dict.get(key).get(f"{skey}_params")
                                for param in score_dict.keys():
                                    if param == "cluster_key" and score_dict.get(param) is not None:
                                        python_params.append("-cluster_cs")
                                        python_params.append(score_dict.get(param))
                            case "co_occurrence":
                                python_params.append("-co_occurrence")
                                score_dict = param_dict.get(key).get(f"{skey}_params")
                                for param in score_dict.keys():
                                    if param == "cluster_key" and score_dict.get(param) is not None:
                                        python_params.append("-cluster_co")
                                        python_params.append(score_dict.get(param))
                                    if param == "interval" and score_dict.get(param) is not None:
                                        python_params.append("-interval")
                                        python_params.append(str(int(score_dict.get(param))))
                                    if param == "n_splits" and score_dict.get(param) is not None:
                                        python_params.append("-n_splits")
                                        python_params.append(str(int(score_dict.get(param))))
                            case "neighborhood_enrichment":
                                python_params.append("-nhood_enrichment")
                                score_dict = param_dict.get(key).get(f"{skey}_params")
                                for param in score_dict.keys():
                                    if param == "cluster_key" and score_dict.get(param) is not None:
                                        python_params.append("-cluster_nhood")
                                        python_params.append(score_dict.get(param))
                                    if param == "library_key" and score_dict.get(param) is not None and score_dict.get(param) != "":
                                        python_params.append("-library_key")
                                        python_params.append(None if score_dict.get(param)=="" else score_dict.get(param))
                                    if param == "n_perms" and score_dict.get(param) is not None:
                                        python_params.append("-n_perms_nhood")
                                        python_params.append(str(int(score_dict.get(param))))
            case "chromVar":
                for ckey in param_dict.get(key).keys():
                    if (not ckey.endswith("_params")) and (param_dict.get(key).get(ckey)):
                        match ckey:
                            case "moranI":
                                multiome_params_py.append("-moranI")
                                score_dict = param_dict.get(key).get(f"{ckey}_params")
                                for param in score_dict.keys():
                                    if param == "n_perms" and score_dict.get(param) is not None:
                                        multiome_params_py.append("-n_perms_autocorr_mI")
                                        multiome_params_py.append(str(int(score_dict.get(param))))
                                    if param == "two_tailed" and score_dict.get(param):
                                        multiome_params_py.append("-two_tailed_mI")
                                    if param == "corr_method" and score_dict.get(param) is not None:
                                        multiome_params_py.append("-corr_method_mI")
                                        multiome_params_py.append(score_dict.get(param))
                            case "gearyC":
                                multiome_params_py.append("-gearyC")
                                score_dict = param_dict.get(key).get(f"{ckey}_params")
                                for param in score_dict.keys():
                                    if param == "n_perms" and score_dict.get(param) is not None:
                                        multiome_params_py.append("-n_perms_autocorr_gC")
                                        multiome_params_py.append(str(int(score_dict.get(param))))
                                    if param == "two_tailed" and score_dict.get(param):
                                        multiome_params_py.append("-two_tailed_gC")
                                    if param == "corr_method" and score_dict.get(param) is not None:
                                        multiome_params_py.append("-corr_method_gC")
                                        multiome_params_py.append(score_dict.get(param))
                            case "differential_motif_activity":
                                multiome_params_py.append("-differential_motif_activity")
                                multiome_params.append("--differential_motif_activity")

            case "liana":
                for lkey in param_dict.get(key).keys():
                    if lkey == "composition_column" and param_dict.get(key).get(lkey) is not None:
                        python_params.append("-cell_comp_key")
                        python_params.append(param_dict.get(key).get(lkey))

    return (python_params, network_params, multiome_params, compute_R_scores, multiome_params_py)

#          job_dir, payload
async def calculate_scores_helper(job_dir, json_dict):
    try:
        # create directory for files created during computations
        base_dir_choices = ["plasmidpoop", "junkDNA420", "kackhaufen1", "dumpase1"]
        base_dir = random.choice(base_dir_choices)

        out_dir = os.path.join(job_dir, base_dir)
        os.makedirs(out_dir, exist_ok=True)

        # create log file
        log_file = os.path.join(out_dir, "calc_scores.log")
        print(log_file)

        if len(json_dict.keys()) == 2 and "params_python_script" in json_dict.keys() and "params_R_script" in json_dict.keys():
            python_params = json_dict.get("params_python_script")
            R_params = json_dict.get("params_R_script")

            with open(log_file, "w") as f:
                f.write("Parameters for python pipeline:" + "\n")
                f.write(str(python_params) + "\n")
                f.write("Parameters for R pipeline:" + "\n")
                f.write(str(R_params) + "\n\n")
        else:
            # get parameters from json_dict
            python_params, R_params, multiome_params, compute_R_scores, multiome_params_py = dict2params(json_dict)

        print(python_params)
        print("these are R_params:")
        print(R_params)

        subprocess.run(["python3", "../backend/calc_python_scores/calc_scores.py",
                            "-outdir", out_dir,
                            "-log", log_file] + python_params,
                            check=True)

        if multiome_params:
            subprocess.run(["Rscript", "../backend/calc_multiome_scores/calc_multiome_scores.R",
                            "--outdir", out_dir,
                            "--log", log_file] + multiome_params,
                            check=True)
            subprocess.run(["python3", "../backend/calc_python_scores/add_to_adata.py",
                            "-indir", out_dir,
                            "-log", log_file,
                            "-multiome"],
                            check=True)

            subprocess.run(["python3", "../backend/calc_multiome_scores/calc_multiome_scores.py",
                            "--dir", out_dir,
                            "--log", log_file] + multiome_params_py,
                            check=True)


        if compute_R_scores:
            subprocess.run(["Rscript", "../backend/calc_R_scores/calc_scores.R",
                            "--dir", out_dir,
                            "--log", log_file] + R_params,
                            check=True)

            subprocess.run(["python3", "../backend/calc_python_scores/add_to_adata.py",
                            "-indir", out_dir,
                            "-log", log_file,
                            "-Rscores"],
                            check=True)

            # delete temporary folders
            for folder in ["expr_info_st", "Rscores_st", "expr_info_tg", "Rscores_tg"]:
                path = os.path.join(out_dir, folder)
                if os.path.exists(path):
                    shutil.rmtree(path)  # removes the whole folder tree


        # finish the log file
        message = f"Finished! Check out the log file and the AnnData object(s) in {out_dir} for details."
        with open(log_file, "a") as f:
            f.write(message + "\n")
        print(message)

        return out_dir

    except subprocess.CalledProcessError as e:
        # error from one of the subprocesses
        err_msg = f"Step failed: {' '.join(e.cmd)} with exit code {e.returncode}"
        with open(log_file, "a") as f:
            f.write(err_msg + "\n")
        print(err_msg)
        return None

    except Exception as e:
        # all other errors
        err_msg = f"Unexpected error: {str(e)}"
        print(err_msg)
        return None


if __name__ == "__main__":

    visium_files = "../backend/datasets_prepro_new"
    for file in os.listdir(visium_files):
        if file.startswith("GSM"):
            file_path = os.path.join(visium_files, file)

            python_params = ["-input", file_path,
                            "-tangram", "-sc_path", "../backend/datasets_prepro_new/Wu_annotated_prepro.h5ad",
                            "-cell_label", "cell_subclass", "-ensembl_col", "ensembl_id", "-feature_col", "feature_type",
                            "-liana", "-cell_comp_key", "celltype_scores",
                            "-moranI", "-gearyC", "-centrality_scores", "-co_occurrence", "-nhood_enrichment",
                            "-R_scores"]
            R_params = ["--tangram",
                        "--sponge_network", "../backend/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv",
                        "--sponge_analysis", "../backend/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv",
                        "--genie_network", "../backend/networks/GENIE3/BRCA/genie3_BRCA_tpm.top_100k.csv",
                        "--aucell", "--gsva", "--ssgsea", "--viper"]

            params_dict = {"params_python_script": python_params,
                        "params_R_script": R_params}

            asyncio.run(calculate_scores_helper(os.path.join("../backend/datasets_scores", file.replace(".h5ad", "")), params_dict))

            # ! needs to be called from frontend directory
