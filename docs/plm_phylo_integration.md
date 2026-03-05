# PLM + Phylogeny integration logic

PhyloFoundry now separates inference into three layers:

1. **Embedding geometry** (`regime_shift`) proposes candidate clade/branch functional regime shifts.
2. **Attention-derived HA sites** (`ha_sites`) nominate residue-level candidates using attention received.
3. **Phylogenetic tests** (ASR + HyPhy + substitution summaries) provide mechanistic testing.

PLM outputs are explicitly treated as hypothesis-generating signals, not direct evidence of selection.

Primary integrated output:
- `summary/site_evidence.tsv`
- `summary/evidence_summary.json`
