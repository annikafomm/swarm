import pandas as pd
from pyscenic.aucell import aucell
from pyscenic.utils import load_motifs
from pyscenic.utils import GeneSignature
import numpy as np
import os
import scanpy as sc
import anndata
import json

def compute_aucell_scores(adata: anndata.AnnData, gene_signatures: list, output_dir: str):
    """
    Compute AUCell scores for gene signatures in the provided AnnData object.

    Parameters:
    - adata: anndata.AnnData object containing gene expression data.
    - gene_signatures: List of GeneSignature objects containing gene sets.
    - output_dir: Directory to save the AUCell scores.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    print(adata.var_names)
    
    adata_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)

    # Compute AUCell scores
    aucell_scores = aucell(adata_df, gene_signatures, num_workers=8)

    print(f"Zeros in AUCell scores: {np.sum(aucell_scores == 0)}")

    
    # Save AUCell scores to AnnData object
    adata.obsm['aucell'] = aucell_scores
    genesets = {str(gs.name): [str(g) for g in list(gs.gene2weight.keys())] for gs in gene_signatures}

    return adata, genesets

def get_signatures_from_genie3(gene_regulatory_csv) -> list:
    """
    Convert a gene regulatory network DataFrame to a list of GeneSignature objects.
    Parameters:
    - gene_regulatory_csv: DataFrame containing the regulatory network with columns ['regulatoryGene', 'targetGene', 'weight'].
    Returns:
    - List of GeneSignature objects.
    """
    if not {'regulatoryGene', 'targetGene', 'weight'}.issubset(gene_regulatory_csv.columns):
        raise ValueError("Input DataFrame must contain 'regulatoryGene', 'targetGene', and 'weight' columns.")
    genesets = []
    for regulator, group in gene_regulatory_csv.groupby("regulatoryGene"):
        gene2weight = dict(zip(group["targetGene"], group["weight"]))
        geneset = GeneSignature(name=regulator, gene2weight=gene2weight)
        genesets.append(geneset)
    return genesets


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute AUCell scores for gene signatures.")
    parser.add_argument("--adata", type=str, required=True, help="Path to the AnnData object.")
    parser.add_argument("--reg_network_csv", type=str, help="Path to the GENIE3/Sponge CSV file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save AUCell scores.")
    args = parser.parse_args()

    # Load AnnData object
    adata = anndata.read_h5ad(args.adata)
    # Load gene regulatory network if provided
    if args.reg_network_csv:
        gene_regulatory_csv = pd.read_csv(args.reg_network_csv)
        # Check if genenames are the same in adata and gene signatures
        adata_genes = set(adata.var_names)
        intersect= set(gene_regulatory_csv["targetGene"]).intersection(adata_genes)
        print(f"Intersecting genes: {len(intersect)} out of {len(set(gene_regulatory_csv['targetGene']))} in the regulatory network and {len(set(adata_genes))} in the AnnData object.")
        gene_signatures = get_signatures_from_genie3(gene_regulatory_csv)
    else:
        raise ValueError("Please provide a path to the GENIE3/Sponge CSV file.")
    

    
    
    # Compute AUCell scores
    adata, gene_sets = compute_aucell_scores(adata, gene_signatures, args.output_dir)
    

    # Save the updated AnnData object with AUCell scores
    file_name = os.path.basename(args.adata).rsplit('.')[0]
   
    adata.write(os.path.join(args.output_dir, f"{file_name}_aucell.h5ad"))
    # Save gene sets to a JSON file
    gene_sets_file = os.path.join(args.output_dir, f"{file_name}_aucell_genesets.json")
    with open(gene_sets_file, 'w') as f:
        
        json.dump(gene_sets, f, indent=4)
    print(f"AUCell scores computed and saved to {args.output_dir}")
    