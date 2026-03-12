#!/usr/bin/env bash

#Usage:
#bash split_faa_by_description.sh input.faa split_proteins

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 input.faa output_dir"
    exit 1
fi

input_fasta="$1"
outdir="$2"

mkdir -p "$outdir"

awk -v outdir="$outdir" '
function trim(s) {
    sub(/^[ \t\r\n]+/, "", s)
    sub(/[ \t\r\n]+$/, "", s)
    return s
}

function sanitize_filename(s,    tmp) {
    s = trim(s)

    # Replace spaces with underscores
    gsub(/[[:space:]]+/, "_", s)

    # Remove parentheses
    gsub(/[()]/, "", s)

    # Replace slashes and backslashes
    gsub(/[\/\\]/, "_", s)

    # Replace other problematic filename characters
    gsub(/[:;,*?"<>|]/, "_", s)

    # Replace plus signs
    gsub(/\+/, "plus", s)

    # Replace ampersands
    gsub(/&/, "and", s)

    # Collapse repeated underscores
    gsub(/_+/, "_", s)

    # Remove leading/trailing underscores or dots
    gsub(/^[_\.]+|[_\.]+$/, "", s)

    if (s == "") s = "unknown_protein"
    return s
}

BEGIN {
    current_file = ""
}

/^>/ {
    header = substr($0, 2)

    # Everything after the first whitespace is treated as the protein description
    desc = header
    sub(/^[^[:space:]]+[[:space:]]*/, "", desc)

    # If there is no description, fall back to the sequence ID
    if (desc == header || desc == "") {
        desc = header
        sub(/[[:space:]].*$/, "", desc)
    }

    fname = sanitize_filename(desc)

    count[fname]++
    if (count[fname] > 1) {
        fname = fname "_" count[fname]
    }

    current_file = outdir "/" fname ".faa"
    print ">" header > current_file
    next
}

{
    if (current_file != "") {
        print >> current_file
    }
}
' "$input_fasta"


