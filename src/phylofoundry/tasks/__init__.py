"""Task modules for the PhyloFoundry pipeline.

Keep package import lightweight: individual task modules are imported lazily
where they are needed to avoid importing optional heavy dependencies at import time.
"""

__all__ = [
    "prep",
    "hmmer",
    "diamond",
    "extract",
    "embed",
    "maape",
    "phylo",
    "post",
    "synteny",
    "codon",
    "hyphy",
    "asr",
    "score_motifs",
    "discover_motifs",
    "curate",
    "aa_composition",
]
