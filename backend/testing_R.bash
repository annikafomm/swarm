#!/bin/bash

R_SCRIPT="./gene_set_ES/calc_scores.R"



IN_DIR="./datasets_scores/GSM6592049_M2_prepro/" 
OUT_DIR="./datasets_scores/GSM6592049_M2_scores/" 

mkdir -p "$OUT_DIR"

echo "Processing $IN_DIR"
Rscript ./gene_set_ES/calc_scores.R --indir ./datasets_scores/GSM6592049_M2_prepro/ --outdir ./datasets_scores/GSM6592049_M2_scores/ --sponge_network ./networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_interactionNetwork.csv --sponge_analysis ./networks/SPONGE/breast_invasive_carcinoma/breast_invasive_carcinoma_networkAnalysis.csv --genie_network ./networks/GENIE3/BRCA/genie3_BRCA_mrn.top_100k.csv --aucell --oe --gsva --ssgsea --viper
echo "Calculated R scores"
echo "Finished! You can find the files with all computed scores here: $OUT_DIR"


powershell -c "(New-Object Media.SoundPlayer 'C:\Windows\Media\Windows Notify.wav').PlaySync()"
echo "All files processed. Check $OUT_DIR/$LOGFILE for details."
# End of script
exit 0