#!/usr/bin/bash

curl -c cookies.txt \
    -X POST http://127.0.0.1:8000/create_session/mopitas && echo ""

curl -b cookies.txt http://127.0.0.1:8000/whoami && echo ""

curl -b cookies.txt \
    -H "Content-Type: application/json" \
    -X POST http://127.0.0.1:8000/read_adata \
    -d '{"path": "../backend/data/adata.h5ad"}' && echo ""

curl -b cookies.txt http://127.0.0.1:8000/obs/cell_type &> /dev/null && echo "obs check successful"

curl -b cookies.txt http://127.0.0.1:8000/var/n_cells &> /dev/null && echo "var check successful"

curl -b cookies.txt http://127.0.0.1:8000/obsm/ligand_receptor_cosine_similarity/LGALS9^PTPRC &> /dev/null && echo "obsm check successful"

curl -b cookies.txt -X POST http://127.0.0.1:8000/delete_session && echo ""
