# First full pipeline QC checklist

Inspect these outputs first:

- `qc/combined/hmmer_bitscore_pre.png`, `hmmer_bitscore_post.png`
- `qc/combined/regime_shift_scores.png`, `regime_shift_pvalues.png`
- `qc/<hmm>/ha_site_counts_per_sequence.png`, `ha_frequency_msa.png`
- `summary/step_status.tsv`, `summary/qc_manifest.tsv`

## Common failure mode: discover finds nothing

1. Verify `summary/ha_sites.tsv` exists and has non-zero `is_ha`.
2. Relax HA thresholds (`ha.call_mode`, `ha.percentile`, or `ha.topk`).
3. Confirm clades are available (`summary/clade_assignment.tsv` / detected clades).
4. Check embedding quality and regime-shift significance.
