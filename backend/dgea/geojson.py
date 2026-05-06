import json

# load data
with open("/workspaces/mopitas-mapra/frontend/public/assets/hexagons.geojson", "r") as f:
    hexagons = json.load(f)

with open("/workspaces/mopitas-mapra/backend/dgea/results/dgea_results.json", "r") as f:
    dgea = json.load(f)

# safety check: ensure hexagons have "meta" field
hexagons.setdefault("meta", {})

# merge DGEA-results in hexagons meta
hexagons["meta"]["dgea"] = dgea["meta"]["dgea"]

# save new file
with open("hexagons_with_dgea.geojson", "w") as f:
    json.dump(hexagons, f, indent=4)

print("Ready: hexagons_with_dgea.geojson created")
