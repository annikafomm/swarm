import scanpy as sc
import mygene
import pandas as pd

adata = sc.read_h5ad("../datasets_prepro/Wu_annotated_prepro.h5ad")

print(adata)


# Initialize mygene client
mg = mygene.MyGeneInfo()

# Get your gene symbols
symbols = list(set(adata.var_names.tolist()))  # or adata.var["gene_symbol"].tolist()

# Query Ensembl
results = mg.querymany(
    symbols,
    scopes="symbol",       # input type
    fields="ensembl.gene,ensembl.type_of_gene,type_of_gene", # output fields
    species="human"        # change to "mouse" if needed
)

# Convert results into dataframe
df = pd.DataFrame(results)
allowed_types = {"miRNA", "protein-coding", "lncRNA"}

def extract_preferred_ids(gene_list: list, allowed_types):
    ids = ""
    gene_type = ""
    for l in gene_list:
        match l.get("type_of_gene"):
            case "miRNA":
                ids = l["gene"]
                gene_type = "miRNA"
                break
            case "protein-coding":
                ids = l["gene"]
                gene_type = "protein_coding"
                break
            case "lncRNA":
                ids = l["gene"]
                gene_type = "lncRNA"
                break

    if len(ids) == 0:
        ids = gene_list[0]["gene"]
        gene_type = gene_list[0]["type_of_gene"]
       
    return ids, gene_type


# Merge with adata.var
mapping = df[["query", "ensembl"]].dropna()
mapping["ensembl_id"] = mapping["ensembl"].apply(
    lambda x: x["gene"] if isinstance(x, dict) else extract_preferred_ids(x, allowed_types)[0] if isinstance(x, list) else None
)
mapping["feature_type"] = mapping["ensembl"].apply(
    lambda x: x["type_of_gene"] if isinstance(x, dict) else extract_preferred_ids(x, allowed_types)[1] if isinstance(x, list) else None
)
mapping = mapping[["query", "ensembl_id", "feature_type"]].drop_duplicates()

# Add to AnnData
adata.var["ensembl_id"] = adata.var_names.map(
    dict(zip(mapping["query"], mapping["ensembl_id"]))
)
adata.var["feature_type"] = adata.var_names.map(
    dict(zip(mapping["query"], mapping["feature_type"]))
)

print(adata)
print(adata.var.head(20))
adata.write_h5ad("../datasets_prepro_new/Wu_annotated_prepro.h5ad")