from typing import List
from pydantic import BaseModel, Field


class TaxonomyItem(BaseModel):
    taxonomy_id: int
    name: str
    family: str


class TaxonomyRequest(BaseModel):
    taxonomy_ids: List[int] = Field(
        min_length=1, description="List of taxonomy ids", examples=[9606, 1, 11292]
    )


class TaxonomyResponse(BaseModel):
    taxonomy: List[TaxonomyItem] = Field(description="List of taxonomy lookup results")

class ClusterRequest(BaseModel):
    sequence_data: dict[str, str] = Field(
        description="Dictionary mapping sequence IDs to raw amino acid sequences",
        examples=[{"seq_01": "MVKV...", "seq_02": "MSKG..."}]
    )
    sequence_identity_threshold: float = Field(
        default=0.5, 
        ge=0.0, 
        le=1.0, 
        description="Sequence identity threshold for clustering (between 0.0 and 1.0)"
    )
    use_linear_clustering: bool = Field(
        default=False, 
        description="If True, uses easy_linclust instead of easy_cluster for massive datasets"
    )


class ClusterResponse(BaseModel):
    clustered_data: dict[str, str] = Field(
        description="Filtered dictionary containing only the representative cluster sequences"
    )