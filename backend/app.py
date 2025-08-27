#from flask import Flask, request, jsonify
import os
import subprocess
import random
import shutil
from datetime import datetime

from backend.calc_python_scores.translate_json import dict2params

#          job_dir, payload
def upload(job_dir, json_dict):
    try:
        # create directory for files created during computations
        base_dir_choices = ["plasmidpoop", "junkDNA420", "kackhaufen1", "dumpase1"]
        base_dir = random.choice(base_dir_choices)

        out_dir = os.path.join(job_dir, base_dir)
        os.makedirs(out_dir, exist_ok=False)

        # create log file
        job_id = json_dict.get('jobId')
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        log_file = os.path.join(out_dir, f"{job_id}_{timestamp}.log")
        print(log_file)

        # get parameters from json_dict
        python_params, R_params = dict2params(json_dict)
        print(python_params)
        print(R_params)

        # Run scripts sequentially
        subprocess.run(["python3", "backend/calc_python_scores/calc_scores.py", 
                        "-outdir", out_dir, 
                        "-log", log_file] + python_params,
                        check=True)

        if R_params:
            subprocess.run(["Rscript", "backend/calc_R_scores/calc_scores.R", 
                            "--dir", out_dir,
                            "--log", log_file] + R_params, 
                            check=True)
        
            subprocess.run(["python3", "backend/calc_python_scores/add_to_adata.py", 
                            "-indir", out_dir,
                            "-log", log_file], 
                            check=True)
        
        # delete temporary folders
        for folder in ["expr_info", "Rscores"]:
            path = os.path.join(out_dir, folder)
            if os.path.exists(path):
                shutil.rmtree(path)  # removes the whole folder tree
                print(f"Deleted {path}")
            else:
                print(f"Skipped {path}, does not exist")

        # finish the log file  
        message = f"Finished! Check out the log file and the AnnData object in {out_dir} for details."
        with open(log_file, "a") as f:
            f.write(message + "\n")
        print(message)

        #return jsonify({"status": "success", "jobId": job_id})
        
    except Exception as e:
        #return jsonify({"status": "error", "message": str(e)}), 500
        print("Error")
