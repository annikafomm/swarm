#!/bin/bash
# filepath: testing_scores.sh


LOGFILE="output.log"
SCRIPT="calc_scores/calc_scores.py"
DIR="datasets_prepro" # "datasets_prepro"

params="_defaults" # _defaults _permsX100 _autocorr_obs_leiden
OUT_DIR="datasets_scores$params" #"outputs/prepro$params"


mkdir -p "$OUT_DIR"

# Empty the log file at the start
> "$OUT_DIR/$LOGFILE"


for file in "$DIR"/G*.h5ad; do
    base=$(basename "$file" .h5ad)

    echo "Processing $file" | tee -a "$OUT_DIR/$LOGFILE"
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment >> "$OUT_DIR/$LOGFILE" 2>&1 ## _defaults
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -n_perms_autocorr 10 -nhood_enrichment -n_perms_nhood 10000 >> "$OUT_DIR/$LOGFILE" 2>&1
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -n_perms_autocorr 10 -attr obs -genes 'leiden' >> "$OUT_DIR/$LOGFILE" 2>&1
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -co_occurrence -n_splits 200 >> "$OUT_DIR/$LOGFILE" 2>&1
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment -filter -normalize >> "$OUT_DIR/$LOGFILE" 2>&1 ## _defaults
    python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment >> "$OUT_DIR/$LOGFILE" 2>&1 ## _defaults
    echo "Finished $OUT_DIR/$base$params.h5ad" | tee -a "$OUT_DIR/$LOGFILE"
    
done

powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $OUT_DIR/$LOGFILE for details." | tee -a "$OUT_DIR/$LOGFILE"
# End of script
exit 0