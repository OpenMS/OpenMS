"""Build target + foreign-like entrapment + decoy FASTAs: every target protein gets one FOREIGN_<id> protein whose
residues are fully shuffled (seeded by the accession), so its tryptic peptides have unrelated lengths and masses, as in
a foreign proteome of equal size. Decoys are full reversals with the DECOY_ prefix, for targets and foreign proteins."""
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
foreign = []
for h, s in targets:
    t = list(s); random.Random(int(hashlib.sha256(h.split()[0].encode()).hexdigest()[:16], 16)).shuffle(t)
    foreign.append(("FOREIGN_" + h, "".join(t)))
with dst.open("w") as f:
    for h, s in targets + foreign: f.write(f">{h}\n{s}\n")
    for h, s in targets + foreign: f.write(f">DECOY_{h}\n{s[::-1]}\n")
print(dst, len(targets), len(foreign), hashlib.sha256(dst.read_bytes()).hexdigest())
