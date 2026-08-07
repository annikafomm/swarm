#!/bin/bash

# Google Drive folder URL (ID: 1j4ngO_QNJONZLfR7SLS-S1No4mZDm-uj)
FOLDER_URL="https://drive.google.com/drive/folders/1j4ngO_QNJONZLfR7SLS-S1No4mZDm-uj"
TEMP_DIR="./temp_downloads"
BACKEND_DIR="./backend/data"
FRONTEND_DIR="./frontend/public/assets"

# Create directories only if they don't exist
mkdir -p "$TEMP_DIR" "$BACKEND_DIR" "$FRONTEND_DIR"

echo "Downloading folder recursively with gdown..."
# test if gdown is installed
if command -v gdown &> /dev/null; then
    gdown --folder "$FOLDER_URL" -O "$TEMP_DIR" --fuzzy
else
    echo "gdown not found. Trying alternative methods..."
    wget -O $TEMP_DIR/adata.h5ad "https://drive.usercontent.google.com/download?id=1WBN1z3zlOBnu58jOvcBeWcerFviHffh9&export=download&confirm=yes"
    wget -O $TEMP_DIR/genie_network_filt.csv "https://drive.usercontent.google.com/download?id=1LqtKT65_jMp0uevbNu3OsKC0bsnQE6yB&export=download&confirm=yes"
    wget -O $TEMP_DIR/hexagons.geojson "https://drive.usercontent.google.com/download?id=1zQyAk-lfY8SdmF4gfvXFcmWBL85dJWWY&export=download&confirm=yes"
    wget -O $TEMP_DIR/sponge_network_smaller.csv "https://drive.usercontent.google.com/download?id=1qHhdBJTTH5kFdAeX6FvwfM4Ji9twes-O&export=download&confirm=yes"
fi

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
