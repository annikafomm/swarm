#!/usr/bin/bash

# Quickly populate the webpage with very stupid random data,
# for debugging purposes only. Overwrites files, so be careful!
# Usage: ./generate_dummy_data.sh GSM6592049_M2.h5ad

BASEDIR=$(dirname $0)

mkdir -p $BASEDIR/backend/data

python $BASEDIR/backend/simulate_scores.py --adata $1 --output $BASEDIR/backend/data

python $BASEDIR/backend/visium_to_geojson.py --adata $BASEDIR/backend/data/adata.h5ad \
    --outpath $BASEDIR/frontend/public/assets/hexagons.geojson
