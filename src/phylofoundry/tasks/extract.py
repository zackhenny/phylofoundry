import os
import glob
from collections import defaultdict
from ..utils.bio import write_fasta

def run_extract(cfg, best_df, fasta_dir, hmm_keep, proteome_seqs, force=False):
    """Build per-HMM FASTA files from competitive best-hit assignments.

    Parameters
    ----------
    cfg : dict
        Pipeline configuration dictionary.
    best_df : pandas.DataFrame
        Competitive best-hit table produced by ``hmmer.best_hits()`` or
        ``diamond.run_diamond()``.  Each row represents the single best HMM
        assignment for a (genome, protein) pair.  When ``phylo.keep_all_hits``
        is ``True`` in the config, the caller should pass the full filtered
        hits table instead so that every hit is retained.
    fasta_dir : str
        Directory where per-HMM ``.faa`` files are written.
    hmm_keep : set or None
        Optional set of HMM names to include.  ``None`` means keep all.
    proteome_seqs : dict
        Nested mapping ``{genome: {protein_id: sequence}}`` used to look up
        amino-acid sequences for each hit.
    force : bool
        When ``True``, overwrite existing per-HMM FASTA files.
    """
    print("\n[extract] Building per-HMM FASTAs from best-hit assignments...")

    hmm_to_seqs = defaultdict(dict)

    if best_df is None or best_df.empty:
        print("WARNING: No hits available to extract sequences.")
        return hmm_to_seqs

    for _, r in best_df.iterrows():
        hmm = r["hmm"]
        if hmm_keep is not None and hmm not in hmm_keep:
            continue

        genome = r["genome"]
        prot = r["protein"]
        seq = proteome_seqs.get(genome, {}).get(prot)
        if not seq:
            continue
        tip = f"{genome}|{prot}"
        hmm_to_seqs[hmm][tip] = seq

    for hmm, seqs in hmm_to_seqs.items():
        out_fp = os.path.join(fasta_dir, f"{hmm}.faa")
        if (not os.path.exists(out_fp)) or force:
            write_fasta(out_fp, seqs)

    return hmm_to_seqs
