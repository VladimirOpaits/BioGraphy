import numpy as np
from core.embedder import BioGraphyEmbedder
from api.api_clients import MyGeneClient, EnsemblClient

def main():
    mg = MyGeneClient()
    ens = EnsemblClient()
    embedder = BioGraphyEmbedder()

    gene_symbol = 'TP53'

    ensembl_id = mg.get_ensembl_id(gene_symbol)
    if not ensembl_id:
        print(f'could not find {gene_symbol}')
        return
    
    dna_sequence = ens.get_sequence_by_id(ensembl_id, chunk_size = 600)
    if not dna_sequence:
        print(f'could not fetch sequence for {ensembl_id}')
        return
    
    vector = embedder.get_vector(dna_sequence)
    print(f"  Vector Snippet: {vector[0][:5]}")

if __name__ == "__main__":
    main()