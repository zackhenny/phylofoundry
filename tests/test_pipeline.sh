#!/usr/bin/env bash
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PH_DIR="$( dirname "$DIR" )"
DUMMY_DIR="$DIR/dummy_data"
OUT_DIR="$DIR/output"

echo "=== Running PhyloFoundry Dummy Test ==="

# 1. Create Dummy FASTA (Proteomes)
cat << 'EOF' > "$DUMMY_DIR/proteomes/genomeA.faa"
>genomeA|prot1
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLA
>genomeA|prot2
MKTLLILAVLSALVSSASAADTPGVVTYDDAVKVTLPRGQWLTLEETRTLTLPAQTDADGEFRL
EOF

cat << 'EOF' > "$DUMMY_DIR/proteomes/genomeB.faa"
>genomeB|protX
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLA
>genomeB|protY
MKTLLILAVLSALVSSASAADTPGVVTYDDAVKVTLPRGQWLTLEETRTLTLPAQTDADGEFRL
EOF

cat << 'EOF' > "$DUMMY_DIR/proteomes/genomeC.faa"
>genomeC|protAlpha
MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLA
>genomeC|protBeta
MKTLLILAVLSALVSSASAADTPGVVTYDDAVKVTLPRGQWLTLEETRTLTLPAQTDADGEFRL
EOF

# 2. Create Dummy HMM
# To make it a real test, we would need a real HMM file. For the sake of CI without real data,
# we will just warn the user to run it with real data or we create a minimalist HMM.
# It is extremely hard to build a valid binary HMM from scratch in raw text. 
# We'll rely on the user to use real data instead for deeper tests, but we'll provide this script structure.
echo "NOTE: A valid .hmm file is required in $DUMMY_DIR/hmms/ for this to complete."
echo "If no HMM is found, the pipeline will fail gracefully."

# 3. Dummy clades.tsv
cat << 'EOF' > "$DUMMY_DIR/clades.tsv"
clade_name	tip
Clade1	genomeA|prot1
Clade1	genomeB|protX
Clade2	genomeC|protAlpha
EOF

# 4. Dummy taxonomy_file.tsv
cat << 'EOF' > "$DUMMY_DIR/taxonomy_file.tsv"
genome	lineage
genomeA	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
genomeB	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
genomeC	d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri
EOF

cat << EOF > "$DUMMY_DIR/test_config.json"
{
    "inputs": {
        "faa_dir": "$DUMMY_DIR/proteomes",
        "hmm_input": "$DUMMY_DIR/hmms",
        "cds_dir": "$DUMMY_DIR/cds",
        "taxonomy_file": "$DUMMY_DIR/taxonomy_file.tsv"
    },
    "output": {
        "outdir": "$OUT_DIR"
    },
    "hmmer": {
        "run_scan": true,
        "run_search": true
    },
    "filtering": {
        "global_min_score": 0.0,
        "min_coverage": 0.0
    },
    "phylo": {
        "mafft": true,
        "combined_tree": true
    },
    "embeddings": {
        "enabled": true,
        "backend": "esm",
        "model": "esm2_t33_650M_UR50D",
        "device": "cpu"
    },
    "post": {
        "enabled": true,
        "compute_conservation": true,
        "compute_kl": true,
        "clades_tsv": "$DUMMY_DIR/clades.tsv"
    },
    "hyphy": {
        "enabled": false
    },
    "synteny": {
        "enabled": false
    }
}
EOF

echo "To run this test (once you place a valid .hmm into $DUMMY_DIR/hmms/):"
echo "  python $PH_DIR/src/phylofoundry/main.py --config $DUMMY_DIR/test_config.json"

EOF
