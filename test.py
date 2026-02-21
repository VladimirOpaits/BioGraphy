import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import requests
import os
import logging

logging.getLogger("transformers").setLevel(logging.ERROR)

class BioGraphyEmbedder:
    def __init__(self, model_name="InstaDeepAI/nucleotide-transformer-500m-human-ref"):
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.tokenizer = AutoTokenizer.from_pretrained(model_name, token = False)
        
        self.model = AutoModel.from_pretrained(
            model_name, 
            dtype=torch.float16,
            trust_remote_code=True,
            token = False
        ).to(self.device)
        self.model.eval()

    def get_vector(self, sequence: str):
        clean_seq = "".join(filter(str.isalpha, sequence.upper()))
        
        inputs = self.tokenizer(
            clean_seq, 
            return_tensors="pt", 
            truncation=True, 
            max_length=1024
        ).to(self.device)
        
        with torch.no_grad():
            outputs = self.model(**inputs)
        
        return outputs.last_hidden_state.mean(dim=1).detach().cpu().numpy()

def get_ensembl_id(gene_symbol):
    print(f"[*] Mapping {gene_symbol} to Ensembl ID...")
    url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=ensembl.gene"
    try:
        res = requests.get(url).json()
        if res.get('hits'):
            return res['hits'][0]['ensembl']['gene']
    except Exception as e:
        print(f"[!] Mapping error: {e}")
    return None

def fetch_dna_sequence(ensembl_id):
    print(f"[*] Fetching sequence for {ensembl_id}...")
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
    except Exception as e:
        print(f"[!] API Fetch error: {e}")
    return None


embedder = BioGraphyEmbedder()

target_gene = "TP53"  # BRCA1, ACTN3, APOE
ens_id = get_ensembl_id(target_gene)

if ens_id:
    dna_seq = fetch_dna_sequence(ens_id)
    
    if dna_seq:
        print(f"[+] Successfully fetched {len(dna_seq)} nucleotides")
        
        fragment = dna_seq[:500] 
        vector = embedder.get_vector(fragment)
        
        print("\n" + "="*60)
        print(f"RESULTS FOR GENE: {target_gene} ({ens_id})")
        print("="*60)
        print(f"Sequence preview: {fragment[:50]}...")
        print(f"Vector shape: {vector.shape}")
        print(f"First 5 vector values: {vector[0][:5]}")
        print("="*60)
    else:
        print("[!] Could not fetch sequence.")
else:
    print("[!] Could not map gene symbol.")