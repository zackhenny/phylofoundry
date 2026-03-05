from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class NormalizedID:
    original: str
    normalized: str
    genome_id: str | None
    protein_id: str


def normalize_sequence_id(raw_id: str, mode: str = "after_last_pipe") -> str:
    sid = str(raw_id or "").strip()
    sid = re.sub(r"\s+", "_", sid)
    if mode == "same":
        return sid
    if mode == "strip_pipe":
        return sid.replace("|", "_")
    if mode == "after_last_pipe":
        return sid.split("|")[-1]
    raise ValueError(f"Unsupported ID mode: {mode}")


def parse_sequence_id(raw_id: str, mode: str = "after_last_pipe") -> NormalizedID:
    sid = str(raw_id or "").strip()
    normalized = normalize_sequence_id(sid, mode=mode)
    parts = sid.split("|")
    genome = parts[0] if len(parts) > 1 else None
    protein = parts[-1]
    return NormalizedID(original=sid, normalized=normalized, genome_id=genome, protein_id=protein)
