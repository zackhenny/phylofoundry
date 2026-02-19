def read_fasta(fp):
    seqs = {}
    h = None
    parts = []
    with open(fp) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(parts)
                h = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if h is not None:
            seqs[h] = "".join(parts)
    return seqs

def write_fasta(fp, records):
    with open(fp, "w") as out:
        for h, s in records.items():
            out.write(f">{h}\n{s}\n")
