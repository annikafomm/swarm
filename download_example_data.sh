#!/bin/bash

# Google Drive folder URL (ID: 1j4ngO_QNJONZLfR7SLS-S1No4mZDm-uj)
FOLDER_URL="https://drive.google.com/drive/folders/1j4ngO_QNJONZLfR7SLS-S1No4mZDm-uj"
TEMP_DIR="./temp_downloads"
BACKEND_DIR="./backend/data"
FRONTEND_DIR="./frontend/public/assets"

# Create directories only if they don't exist
mkdir -p "$TEMP_DIR" "$BACKEND_DIR" "$FRONTEND_DIR"

echo "Downloading folder recursively with gdown..."
# Install gdown if needed: pip install gdown
gdown --folder "$FOLDER_URL" -O "$TEMP_DIR" --fuzzy

echo "Moving files to targets..."
# Move data files to backend/data
find "$TEMP_DIR" -name "adata.h5ad" -exec mv {} "$BACKEND_DIR/" \;
find "$TEMP_DIR" -name "genie_network_filt.csv" -exec mv {} "$BACKEND_DIR/" \;
find "$TEMP_DIR" -name "sponge_network_smaller.csv" -exec mv {} "$BACKEND_DIR/" \;

# Move hexagons.geojson to frontend/public/assets
find "$TEMP_DIR" -name "hexagons.geojson" -exec mv {} "$FRONTEND_DIR/" \;

# Optional: Clean up temp directory (uncomment if desired)
# rm -rf "$TEMP_DIR"

echo "Done! Check backend/data and frontend/public/assets."
