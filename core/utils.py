def get_local_index(global_pos, gene_start, gene_end, strand):
    if strand == 1:
        return global_pos - gene_start
    else:
        return gene_end - global_pos
    
def validate_nucleotide_base(v: str) -> str:
    v_upper = v.upper()
    if v_upper not in 'ACGT':
        raise ValueError(f"Invalid DNA base: {v}")
    return v_upper

def validate_strand(cls, v: int) -> int:
        if v not in (1, -1):
            raise ValueError("Strand must be 1 or -1")
        return v