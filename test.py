import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np

class BioGraphyEmbedder:
    def __init__(self, model_name="InstaDeepAI/nucleotide-transformer-500m-human-ref"):
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(
            model_name, 
            torch_dtype=torch.float16,
            trust_remote_code=True
        ).to(self.device)
        self.model.eval()

    def get_vector(self, sequence: str):
        inputs = self.tokenizer(sequence.upper(), return_tensors="pt", truncation=True, max_length=1024).to(self.device)
        with torch.no_grad():
            outputs = self.model(**inputs)
        return outputs.last_hidden_state.mean(dim=1).detach().cpu().numpy()

embedder = BioGraphyEmbedder()

original = "GCTGTTGCT" * 10 

silent_seq = "GCCGTTGCT" + ("GCTGTTGCT" * 9)

charge_swap_seq = "GCTGATGCT" + ("GCTGTTGCT" * 9)

nonsense_seq = "TAAGTTGCT" + ("GCTGTTGCT" * 9)


def run_test(name, seq1, seq2):
    v1 = embedder.get_vector(seq1)
    v2 = embedder.get_vector(seq2)
    sim = cosine_similarity(v1, v2)[0][0]
    dist = 1 - sim
    print(f"{name:25} | Similarity: {sim:.6f} | Distance: {dist:.6f}")

print("\n" + "="*60)
print("BIOGRAPHY GENOMIC SENSITIVITY REPORT (RTX 4050)")
print("="*60)

run_test("1. Silent (No Change)", original, silent_seq)
run_test("2. Charge Swap (+ to -)", original, charge_swap_seq)
run_test("3. Nonsense (Stop)", original, nonsense_seq)

print("="*60)