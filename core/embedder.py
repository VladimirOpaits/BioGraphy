import os
import logging
import torch
from transformers import AutoTokenizer, AutoModel
from scipy.spatial.distance import cosine  # type: ignore
from .models import MutationPair

os.environ['TRANSFORMERS_VERBOSITY'] = 'error'
logging.getLogger('transformers').setLevel(logging.ERROR)

class BioGraphyEmbedder:
    def __init__(self, model_name = "InstaDeepAI/nucleotide-transformer-500m-human-ref"):
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'

        self.tokenizer = AutoTokenizer.from_pretrained(model_name, token = False)
        self.model = AutoModel.from_pretrained(
            model_name, 
            dtype=torch.float16, 
            trust_remote_code = True, 
            token = False
            ).to(self.device)
        self.model.eval()

    def get_vector(self, sequence: str):
        clean_seq = ''.join(filter(str.isalpha, sequence.upper()))
        inputs = self.tokenizer(
            clean_seq,
            return_tensors = 'pt',
            padding=True,
            truncation = True,
            max_length = self.tokenizer.model_max_length
        ).to(self.device)

        with torch.no_grad():
            outputs = self.model(**inputs)

        embeddings = outputs.last_hidden_state

        embeddings = embeddings.to(torch.float32)
        mask = inputs['attention_mask'].unsqueeze(-1).expand(embeddings.size()).to(torch.float32)

        sum_embeddings = torch.sum(embeddings * mask, 1)
        sum_mask = torch.clamp(mask.sum(1), min=1e-9)
        mean_pooled = sum_embeddings / sum_mask

        return mean_pooled.detach().cpu().numpy().flatten()


    def get_mutation_distance(self, pair: MutationPair) -> float:
        v_ref = self.get_vector(pair.ref)
        v_alt = self.get_vector(pair.alt)

        return float(cosine(v_ref.flatten(), v_alt.flatten()))
