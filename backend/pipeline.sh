#!/bin/bash
# filepath: backend/testing_scores.sh
# Purpose: Run Python and R score calculations on a given dataset and log results

set -euo pipefail  # Stop on errors, undefined variables, and pipeline errors

# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input paths and parameter values
UPLOAD_DIR="results"
INPUT="$SCRIPT_DIR/datasets_prepro/GSM6592049_M2_prepro.h5ad"

SPONGE_NETWORK="$SCRIPT_DIR/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv"
SPONGE_ANALYSIS="$SCRIPT_DIR/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv"

# Extract base name of input file
BASE_NAME=$(basename "$INPUT" .h5ad)

# Random folder name from predefined choices
BASE_DIR_CHOICES=("plasmidpoop" "junkDNA420" "kackhaufen1" "dumpase1")
BASE_DIR=${BASE_DIR_CHOICES[$RANDOM % ${#BASE_DIR_CHOICES[@]}]}

# Create output directory
OUT_DIR="$SCRIPT_DIR/$UPLOAD_DIR/$BASE_DIR"
mkdir -p "$OUT_DIR"

# Log file
LOGFILE="$OUT_DIR/scores.log"
> "$LOGFILE"

# Run Python score calculation
python "$SCRIPT_DIR/calc_python_scores/calc_scores.py" \
    -input "$INPUT" \
    -outdir "$OUT_DIR" \
    -log "$LOGFILE" \
    -moranI

# Run R score calculation
Rscript "$SCRIPT_DIR/calc_R_scores/calc_scores.R" \
    --dir "$OUT_DIR" \
    --log "$LOGFILE" \
    --sponge_network "$SPONGE_NETWORK" \
    --sponge_analysis "$SPONGE_ANALYSIS" \
    --ensembl_col ensemble_id \
    --aucell

# Add R scores to AnnData
python "$SCRIPT_DIR/calc_scores/add_to_adata.py" \
    -indir "$OUT_DIR" \
    -log "$LOGFILE"

# Optional cleanup
# rm -rf "$OUT_DIR/expr_info"
# rm -rf "$OUT_DIR/Rscores"

echo "Finished! Check out the log file and the AnnData object in $OUT_DIR for details." | tee -a "$LOGFILE"

exit 0
