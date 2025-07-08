#!/bin/bash


LOGFILE="output.log"
SCRIPT="preprocessing/preprocessing.py"
mt_ensembl_file="preprocessing/mt_ensembl_ids.txt"
DIR="datasets"

OUT_DIR="datasets_prepro_igraph"


mkdir -p "$OUT_DIR"

# Empty the log file at the start
> "$OUT_DIR/$LOGFILE"


for file in "$DIR"/*.h5ad; do
    base=$(basename "$file" .h5ad)

    echo "Processing $file" | tee -a "$OUT_DIR/$LOGFILE"
    python "$SCRIPT" -input "$file" -output "$OUT_DIR/$base-prepro.h5ad" -mt_ensembl_ids "$mt_ensembl_file" >> "$OUT_DIR/$LOGFILE" 2>&1
    echo "Finished $OUT_DIR/$base-prepro.h5ad" | tee -a "$OUT_DIR/$LOGFILE"
    
done

powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $OUT_DIR/$LOGFILE for details." | tee -a "$OUT_DIR/$LOGFILE"
# End of script
exit 0