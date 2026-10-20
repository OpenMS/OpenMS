"""Write the frozen #10364 HYE target-decoy FASTA: targets (header kept), then
one full-sequence reversal per target with DECOY_ prefix (prepare_databases.py).
The result must have sha256 efbf5ac224bf51fca252f58cd55df5e2aad7eddd824d0dee9959b086e147bb59."""
import hashlib, sys


def read_fasta(path):
    entries, head, seq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if head is not None:
                    entries.append((head, "".join(seq)))
                head, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if head is not None:
        entries.append((head, "".join(seq)))
    return entries


def split_header(h):
    parts = h.split(None, 1)
    return parts[0], (parts[1] if len(parts) > 1 else "")


if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    entries = [split_header(h) + (s,) for h, s in read_fasta(src)]
    with open(out, "w") as f:
        for ident, desc, s in entries:
            f.write(">" + ident + " " + desc + "\n" + s + "\n")  # as prepare_databases.py (trailing blank kept)
        for ident, desc, s in entries:
            f.write(">DECOY_" + ident + "\n" + s[::-1] + "\n")
    print(len(entries), "targets", hashlib.sha256(open(out, "rb").read()).hexdigest())
