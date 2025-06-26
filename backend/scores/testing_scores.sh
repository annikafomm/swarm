#!/bin/bash
# filepath: testing_scores.sh


LOGFILE="output.log"
SCRIPT="calc_scores/calc_scores.py"
DIR="datasets_prepro"

params="_permsX10" # _defaults _permsX100
OUT_DIR="outputs/prepro$params"


mkdir -p "$OUT_DIR"

# Empty the log file at the start
> "$OUT_DIR/$LOGFILE"


for file in "$DIR"/*.h5ad; do
    base=$(basename "$file" .h5ad)

    echo "Processing $file" | tee -a "$OUT_DIR/$LOGFILE"
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment >> "$OUT_DIR/$LOGFILE" 2>&1 ## _defaults
    python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -n_perms_autocorr 10 >> "$OUT_DIR/$LOGFILE" 2>&1
    #python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base$params.h5ad" -moranI -gearyC -genes 'ENSG00000197063' 'ENSG00000247077' 'ENSG00000121274' 'ENSG00000197106' 'ENSG00000101096' 'ENSG00000161277' 'ENSG00000110921' 'ENSG00000149798' 'ENSG00000254676' 'ENSG00000165905' 'ENSG00000162572' 'ENSG00000223546' 'ENSG00000275052' 'ENSG00000143149' 'ENSG00000269938' 'ENSG00000141294' 'ENSG00000127955' 'ENSG00000168925' 'ENSG00000118495' 'ENSG00000188368' >> "$OUT_DIR/$LOGFILE" 2>&1
    echo "Finished $OUT_DIR/$base$params.h5ad" | tee -a "$OUT_DIR/$LOGFILE"
done

powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $OUT_DIR/$LOGFILE for details." | tee -a "$OUT_DIR/$LOGFILE"
# End of script
exit 0