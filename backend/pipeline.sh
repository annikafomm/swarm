#!/bin/bash
# filepath: backend/testing_scores.sh

# Get the directory of the current script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT="$SCRIPT_DIR/datasets_prepro/GSM6592049_M2_prepro.h5ad" 
sponge_network="$SCRIPT_DIR/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv"
sponge_analysis="$SCRIPT_DIR/networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv"

base=$(basename "$INPUT" .h5ad)

params="defaults_small" 
OUT_DIR="$SCRIPT_DIR/results/${base}_$params" 

mkdir -p "$OUT_DIR"

LOGFILE_PATH="$OUT_DIR/scores_pipeline.log"
# Empty the log file at the start
> "$LOGFILE_PATH"

echo "Processing $INPUT" | tee -a "$LOGFILE_PATH"

cd "$SCRIPT_DIR/calc_scores"
python calc_scores.py -input "$INPUT" -outdir "$OUT_DIR" -moranI >> "$LOGFILE_PATH" 2>&1 
cd "$SCRIPT_DIR"

echo "Calculated squidpy scores" | tee -a "$LOGFILE_PATH"

# Change to the R project directory
cd "$SCRIPT_DIR/gene_set_ES"
Rscript calc_scores.R --dir "$OUT_DIR" --log "$LOGFILE_PATH" --sponge_network "$sponge_network" --sponge_analysis "$sponge_analysis" --ensembl_col ensemble_id --aucell
cd "$SCRIPT_DIR"

echo "Calculated R scores" | tee -a "$LOGFILE_PATH"

python "$SCRIPT_DIR/calc_scores/add_to_adata.py" -indir "$OUT_DIR"

echo "Finished! You can find the file with all computed scores here: $OUT_DIR/${base}_scores.h5ad" | tee -a "$LOGFILE_PATH"


powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $LOGFILE_PATH for details." | tee -a "$LOGFILE_PATH"
# End of script
exit 0