import mygene
import pandas as pd

def flatten(entry):
    # Flatten the nested 'genomic_pos' dictionary for the DataFrame row
    flat_entry = entry.copy()
    for k, v in entry['genomic_pos'].items():
        flat_entry[f'genomic_pos_{k}'] = v
    del flat_entry['genomic_pos']

    df = pd.DataFrame([flat_entry])
    return df

mg = mygene.MyGeneInfo()

if True:
    # Query all human mitochondrial genes by chromosome
    results = mg.query('symbol:MT-*', species='human', fields='ensembl.gene, symbol', size=1000)
    print("All matches are returned: ", len(results['hits']) < 1000)

    # write the ensembl IDs to a set
    mt_ensembl_ids = set()
    for hit in results['hits']:
        ensg = hit.get('ensembl')
        if isinstance(ensg, list):
            mt_ensembl_ids.update(e['gene'] for e in ensg)
        elif isinstance(ensg, dict):
            mt_ensembl_ids.add(ensg['gene'])


    with open("mt_ensembl_ids.txt", "w") as f:
        for item in mt_ensembl_ids:
            f.write(f"{item}\n")
else:
    df = pd.read_csv(r"C:\RebEll\Computer\BioinformatikMaster\Semester8\MaPra\mopitas-mapra\backend\scores\preprocessing\xenium_var.csv", index_col=0)

    symbols = df.index.tolist()  # your list of gene symbols
    ensembl_ids = df['gene_ids'].tolist()
    results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields=['symbol', 'genomic_pos'], species='human')
    print(len(results))
    print(flatten(results[0]))
    #mito_genes = [r['query'] for r in results if r.get('chromosome') == 'MT']
    #print(mito_genes)
    
