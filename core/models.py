from pydantic import BaseModel, Field, field_validator
from typing import Tuple, List
from .utils import validate_nucleotide_base, validate_strand

class GeneData(BaseModel):
    ensembl_id: str
    symbol: str = Field(..., description="Common gene name, e.g., TP53")
    biotype: str = Field(default="protein_coding")
    sequence: str
    strand: int
    start: int
    end: int
    exons: List[Tuple[int, int]] = Field(default=[], description="List of exon coordinates")

    @field_validator('strand')
    @classmethod
    def check_strand(cls, v: int) -> int:
        return validate_strand(cls, v)

class MutationPair(BaseModel):
    ref: str = Field(..., description="Reference sequence window")  # ... for required fields
    alt: str = Field(..., description="Alternative (mutated) sequence window")   #using str instead of seq for robustness
    original_base: str = Field(..., min_length=1, max_length=1)
    new_base: str = Field(..., min_length=1, max_length=1)
    relative_idx: int
    window_range: Tuple[int, int]

    @field_validator('new_base', 'original_base')
    @classmethod  #validator must be class method
    def check_base(cls, v: str) -> str:
        return validate_nucleotide_base(v)