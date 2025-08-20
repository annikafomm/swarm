#!/bin/bash
# filepath: testing_scores.sh


LOGFILE="output.log"
PY_SCRIPT="./calc_scores/calc_scores.py"
#PY_SCRIPT2="./calc_scores/calc_scores.py"
R_SCRIPT="./gene_set_ES/calc_scores.R"



INPUT="../datasets_prepro/GSM6592049_M2_prepro.h5ad" 

params="defaults" 
OUT_DIR="datasets_scores/$params" 

mkdir -p "$OUT_DIR"

# Empty the log file at the start
> "$OUT_DIR/$LOGFILE"


base=$(basename "$INPUT" .h5ad)
echo "$base"

echo "Processing $INPUT" | tee -a "$OUT_DIR/$LOGFILE"
python "$PY_SCRIPT" -input "$INPUT" -output "$OUT_DIR/$base-$params.h5ad" -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment >> "$OUT_DIR/$LOGFILE" 2>&1 
echo "Calculated squidpy scores" | tee -a "$OUT_DIR/$LOGFILE"
#python "$R_SCRIPT" ...
echo "Calculated R scores" | tee -a "$OUT_DIR/$LOGFILE"
#python "$PY_SCRIPT2" ...
echo "Finished! You can find the file with all computed scores here: $OUT_DIR/$base-$params.h5ad" | tee -a "$OUT_DIR/$LOGFILE"


powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $OUT_DIR/$LOGFILE for details." | tee -a "$OUT_DIR/$LOGFILE"
# End of script
exit 0