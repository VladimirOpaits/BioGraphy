from Bio.Seq import Seq
from core.models import MutationPair

class MutationProcessor:
    def __init__(self, window_size=512):
        self.window_size = window_size

    def create_mutation_pair(self, sequence, local_index, new_base):
        full_seq = Seq(sequence)
        seq_len = len(full_seq)
        
        if self.window_size > seq_len:
            start, end = 0, seq_len
        else:
            start = local_index - self.window_size // 2
            end = start + self.window_size
            
            if start < 0:
                start, end = 0, self.window_size
            elif end > seq_len:
                end = seq_len
                start = seq_len - self.window_size

        ref_window = full_seq[start:end]
        rel_idx = local_index - start
        
        original_base = str(ref_window[rel_idx])
        if original_base == new_base:
            return None
            
        alt_list = list(ref_window)
        alt_list[rel_idx] = new_base
        alt_window = Seq("".join(alt_list))
        
        return MutationPair(
            ref=str(ref_window),
            alt=str(alt_window),
            original_base=original_base,
            new_base=new_base,
            relative_idx=rel_idx,
            window_range=(start, end)
        )

    def get_protein_impact(self, ref_seq, alt_seq):
        ref_prot = Seq(ref_seq).translate(stop_symbol="*", to_stop=False)
        alt_prot = Seq(alt_seq).translate(stop_symbol="*", to_stop=False)
        return str(ref_prot), str(alt_prot)