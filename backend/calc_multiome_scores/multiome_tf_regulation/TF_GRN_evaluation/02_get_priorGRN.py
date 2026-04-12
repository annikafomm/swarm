import pandas as pd
import requests
from io import StringIO
import argparse

parser = argparse.ArgumentParser(description="Generate prior GRN from GO genes")
parser.add_argument("go_genes_path", help="Path to gene module CSV file")
parser.add_argument("output_path", help="Path to output priorGRN CSV file")
args = parser.parse_args()

go_genes_path = args.go_genes_path
output_path = args.output_path

go_genes = pd.read_csv(go_genes_path)["SYMBOL"].tolist()

url = "https://omnipathdb.org/interactions/"
params = {
    "genesymbols": 1,
    "datasets": "collectri",
    "organisms": 9606,
    "targets": ",".join(go_genes),
    "fields": "sources,references"
}

r = requests.get(url, params=params, timeout=60)
r.raise_for_status()

raw = pd.read_csv(StringIO(r.text), sep="\t")

priorGRN = (
    raw[["target_genesymbol", "source_genesymbol"]]
    .rename(columns={"target_genesymbol": "gene", "source_genesymbol": "TF"})
    .dropna()
    .drop_duplicates()
    .sort_values(["gene", "TF"])
    .reset_index(drop=True)
)

# save priorGRN
priorGRN.to_csv(output_path, index=False)

