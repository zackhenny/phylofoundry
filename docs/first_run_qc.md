# First-run QC checklist

After running the pipeline, inspect:

1. `summary/step_status.tsv` for failed stages.
2. `summary/qc_manifest.tsv` for per-HMM hit counts and basic metrics.
3. `qc/combined/hits_bitscore.png` and `qc/combined/hits_qcov.png`.
4. `summary/id_audit.tsv` for ID normalization consistency.
5. `discover/discover_profile_failures.tsv` if discovery failed due to low profile count.

If motif discovery returns nothing, check logs in `logs/phylofoundry.log` and lower discovery thresholds or validate ESM runtime/device settings.
