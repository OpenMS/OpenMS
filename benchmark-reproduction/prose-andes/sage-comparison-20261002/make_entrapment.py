"""Build target + shuffled-entrapment + decoy FASTAs (r = 1) from a benchmark target-decoy FASTA.

Every target protein gets one entrapment protein ENTRAP_<id>: K, R and P stay at their positions and the other residues
are shuffled within each segment between two of them, so cleavage sites, peptide lengths and peptide masses are
preserved (paired-mass entrapment, as in FDRBench). Decoys are full reversals with the DECOY_ prefix, the method of the
benchmark databases, for targets and entrapment alike. Deterministic: the shuffle is seeded by the protein accession."""
import hashlib, random, sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
targets, header, seq = [], None, []
for line in src.read_text().splitlines():
    if line.startswith(">"):
        if header and not header.startswith("DECOY_"): targets.append((header, "".join(seq)))
        header, seq = line[1:], []
    else:
        seq.append(line.strip())
if header and not header.startswith("DECOY_"): targets.append((header, "".join(seq)))


def shuffled(acc, s):
    rng = random.Random(int(hashlib.sha256(acc.encode()).hexdigest()[:16], 16))
    out, seg = [], []
    for ch in s + "K":  # sentinel flushes the last segment
        if ch in "KRP":
            for _ in range(10):
                t = seg[:]; rng.shuffle(t)
                if t != seg or len(set(seg)) < 2: break
            out += t; out.append(ch); seg = []
        else:
            seg.append(ch)
    return "".join(out[:-1])


with dst.open("w") as f:
    entrap = [("ENTRAP_" + h, shuffled(h.split()[0], s)) for h, s in targets]
    for h, s in targets + entrap: f.write(f">{h}\n{s}\n")
    for h, s in targets + entrap: f.write(f">DECOY_{h}\n{s[::-1]}\n")
print(dst, "targets", len(targets), "entrapment", len(entrap), "sha256", hashlib.sha256(dst.read_bytes()).hexdigest())
