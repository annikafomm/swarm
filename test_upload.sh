#!/bin/bash
set -euo pipefail

SPATIAL_FILE="backend/data/adata.h5ad"
MULTIOME_FILE="backend/data/multiome_data/healthy_breast_preprocessed_RNA_peaks_only_7_8.rds"
FRAGMENTS_FILE="backend/data/multiome_data/GSM7821196_KTBpool6_Ashkenazi_jew_atac_fragments.tsv.gz"
FRAGMENTS_TBI="backend/data/multiome_data/GSM7821196_KTBpool6_Ashkenazi_jew_atac_fragments.tsv.gz.tbi"

BASE="http://localhost:3005"
COOKIE_JAR="cookies.txt"

# 1) Create session + store cookie
curl -sS -c "$COOKIE_JAR" -X POST "$BASE/create_session/testuser" >/dev/null

# (optional) check whoami using that cookie
curl -sS -b "$COOKIE_JAR" "$BASE/whoami"; echo

# 2) Upload with cookie + matching field names
curl -X POST "$BASE/api/upload" \
  -b "$COOKIE_JAR" \
  -F "email=test@example.com" \
  -F "dataset=Visium" \
  -F "spatial_h5ad=@$SPATIAL_FILE" \
  -F "spatial_normalization=false" \
  -F "spatial_filtering=false" \
  -F "use_tangram=false" \
  -F "use_multiome=true" \
  -F "multiome_rds=@$MULTIOME_FILE" \
  -F "fragments_tsv_gz=@$FRAGMENTS_FILE" \
  -F "fragments_tsv_gz_tbi=@$FRAGMENTS_TBI" \
  -F "score_network=false" \
  -F "score_squidpy=false" \
  -F "score_liana_plus=false" \
  -F "score_chromVar=false" \
  -F "score_differential_motif_activity=false" \
  -F "score_motif_enrichment=false" \
  -F "score_FootprintingBias=true" \
  -F "genome=hg38" \
  --progress-bar
echo
