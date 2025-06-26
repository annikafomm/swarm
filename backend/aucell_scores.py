import pandas as pd
from pyscenic.aucell import aucell
from pyscenic.utils import load_motifs
from pyscenic.utils import GeneSignature
import numpy as np
import os
import scanpy as sc
import anndata

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
    

    # Compute AUCell scores
    aucell_scores = aucell(adata.X, gene_signatures, num_workers=8)

    
    # Save AUCell scores to AnnData object
    adata.obsm['aucell'] = aucell_scores
    anndata.uns["genesets"] = {gs.name: list(gs.gene2weight.keys()) for gs in gene_signatures}
    # Save the updated AnnData object
    adata.write(os.path.join(output_dir, 'anndata_with_aucell_scores.h5ad'))

def get_signatures_from_genie3(gene_regulatory_csv) -> list:
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
        gene_signatures = get_signatures_from_genie3(gene_regulatory_csv)
    else:
        raise ValueError("Please provide a path to the GENIE3/Sponge CSV file.")
    
    # Compute AUCell scores
    compute_aucell_scores(adata, gene_signatures, args.output_dir)
    print(f"AUCell scores computed and saved to {args.output_dir}")
    