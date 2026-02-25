import numpy as np
from core.embedder import BioGraphyEmbedder
from api.api_clients import MyGeneClient, EnsemblClient
from sklearn.metrics.pairwise import cosine_similarity

def run_mutation_test(name, original_seq, mutant_seq, embedder):
    v_orig = embedder.get_vector(original_seq)
    v_mutant = embedder.get_vector(mutant_seq)
    
    sim = cosine_similarity(v_orig, v_mutant)[0][0]
    dist = 1 - sim
    
    print(f"Test: {name}")
    print(f"  Similarity: {sim:.6f}")
    print(f"  Distance  : {dist:.6f}")
    print("-" * 30)


mg = MyGeneClient()
ens = EnsemblClient()
embedder = BioGraphyEmbedder()

gene_symbol = 'TP53'
ensembl_id = mg.get_ensembl_id(gene_symbol)
    
original_dna = ens.get_sequence_by_id(ensembl_id, chunk_size=600)
    
print(f"\nBIOGRAPHY SENSITIVITY TEST: {gene_symbol}")
print("="*40)
    
if original_dna:
    list_dna = list(original_dna)
    list_dna[300] = 'A' if list_dna[300] != 'A' else 'C'
    snv_dna = "".join(list_dna)
    run_mutation_test("Single Point Mutation", original_dna, snv_dna, embedder)

    deletion_dna = original_dna[:300] + original_dna[320:]
    run_mutation_test("20bp Deletion (Frameshift)", original_dna, deletion_dna, embedder)
