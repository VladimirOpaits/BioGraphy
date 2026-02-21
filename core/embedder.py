import os
import logging
import torch
from transformers import AutoTokenizer, AutoModel

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
        input = self.tokenizer(
            clean_seq,
            return_tensors = 'pt',
            truncation = True,
            max_length = 1024
        ).to(self.device)

        with torch.no_grad():
            outputs = self.model(**input)

        return outputs.last_hidden_state.mean(dim = 1).detach().cpu().numpy()
