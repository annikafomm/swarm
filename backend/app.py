#from flask import Flask, request, jsonify
import os
import subprocess
import random
import shutil
from datetime import datetime

from calc_python_scores.translate_json import dict2params

#          job_dir, payload
def upload(job_dir, json_dict):
    try:
        # create directory for files created during computations
        base_dir_choices = ["plasmidpoop", "junkDNA420", "kackhaufen1", "dumpase1"]
        base_dir = random.choice(base_dir_choices)

        out_dir = os.path.join(job_dir, base_dir)
        os.makedirs(out_dir, exist_ok=True)

        # create log file
        job_id = json_dict.get('jobId')
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file = os.path.join(out_dir, f"{job_id}_{timestamp}.log")
        print(log_file)
        
        # get parameters from json_dict
        #python_params, R_params = dict2params(json_dict)
        python_params = ["-input", "./datasets_prepro_new/GSM6592049_M2_prepro.h5ad", 
                         "-tangram", "-sc_path", "./datasets_prepro/Wu_annotated_prepro.h5ad",
                         "-liana", "-cell_comp_key", "celltype_scores",
                         "-moranI", "-gearyC", "-centrality_scores", "-co_occurrence", "-nhood_enrichment",
                         "-R_scores"]
        R_params = ["--tangram", 
                    "--sponge_network", "./networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv", "--sponge_analysis", "./networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv", "--ensembl_col", "ensemble_id", 
                    "--genie_network", "./networks/GENIE3/BRCA/genie3_BRCA_tpm.top_100k.csv", 
                    "--aucell", "--gsva", "--ssgsea", "--viper"]
        print(python_params)
        print(R_params)
        
        # Run scripts sequentially
        subprocess.run(["python3", "calc_python_scores/calc_scores.py",
                        "-outdir", out_dir,
                        "-log", log_file] + python_params,
                        check=True)
        
        if R_params:
            subprocess.run(["Rscript", "calc_R_scores/calc_scores.R",
                            "--dir", out_dir,
                            "--log", log_file] + R_params,
                            check=True)

            subprocess.run(["python3", "calc_python_scores/add_to_adata.py",
                            "-indir", out_dir,
                            "-log", log_file],
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

        return {"status": "success", "jobId": job_id}

    except subprocess.CalledProcessError as e:
        # error from one of the subprocesses
        err_msg = f"Step failed: {' '.join(e.cmd)} with exit code {e.returncode}"
        with open(log_file, "a") as f:
            f.write(err_msg + "\n")
        print(err_msg)
        return {"status": "error", "message": err_msg}

    except Exception as e:
        # all other errors
        err_msg = f"Unexpected error: {str(e)}"
        print(err_msg)
        return {"status": "error", "message": err_msg}

if __name__ == "__main__":
    upload("./uploads/job_0001", {'jobId':'job_0001'})
