#!/usr/bin/env bash
# Downloads and checks the input FASTAs of build_entrapment_dbs.sh.
#
# usage: prepare_inputs.sh <work dir>
#
# Writes <work dir>/data/{hy,hye,human}_targets.fasta: the target part of the #10364 databases, cut from the ProteoBench
# HYE FASTA exactly as the frozen prepare_databases.py does (its target-decoy files are checked against databases.json),
# and <work dir>/data/arath_sp.fasta (reviewed Arabidopsis thaliana, UP000006548).
set -euo pipefail
D=$1/data
mkdir -p "$D"

curl -fsS -o "$D/hye_source.fasta" \
  https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/ProteoBenchFASTA_MixedSpecies_HYE.fasta
echo "d9ac434d88492c10c8e9a587ee7dbc9480fa0995fa07a6ba35a7da8abf39aa25  $D/hye_source.fasta" | sha256sum -c -

python3 -I - "$D" <<'EOF'
import hashlib, sys
from pathlib import Path
data = Path(sys.argv[1])
entries, hdr, seq = [], None, []
for line in (data / "hye_source.fasta").read_text().splitlines():
    if line.startswith(">"):
        if hdr is not None: entries.append((hdr, "".join(seq)))
        hdr, seq = line[1:], []
    else: seq.append(line.strip())
entries.append((hdr, "".join(seq)))
def split(h):
    h = h.strip(); parts = h.split(None, 1)
    return parts[0], (parts[1] if len(parts) > 1 else "")
E = [(*split(h), s) for h, s in entries]
frozen = {"hy": "35614831eddb913fab678c63461d6b4cb44ca73f5aebe2ead1a00020d791def6",
          "hye": "efbf5ac224bf51fca252f58cd55df5e2aad7eddd824d0dee9959b086e147bb59",
          "human": "6f86369103a9003560661f559cbd5bd1c37e55187670f662d69f91f81359a3ab"}
for db in ["hy", "hye", "human"]:
    es = E
    if db == "hy": es = [e for e in E if e[0].endswith(("_HUMAN", "_YEAST")) or "Cont_" in e[0]]
    if db == "human": es = [e for e in E if e[0].endswith("_HUMAN") or "Cont_" in e[0]]
    targets = "".join(">" + i + " " + d + "\n" + s + "\n" for i, d, s in es)
    td = targets + "".join(">DECOY_" + i + "\n" + s[::-1] + "\n" for i, d, s in es)
    if hashlib.sha256(td.encode()).hexdigest() != frozen[db]:
        sys.exit(f"{db}: target-decoy database differs from the frozen #10364 one")
    (data / f"{db}_targets.fasta").write_text(targets)
    print(db, len(es), "targets; matches the frozen #10364 database")
EOF

curl -fsS -o "$D/arath_sp.fasta" \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000006548%29+AND+%28reviewed%3Atrue%29"
# UniProt release 2026_03, as in #10629; a later release gives another file
echo "e20271b5302eaa14a07a21f05910e16a15ea3c80f49e242bfe14774627a38e8e  $D/arath_sp.fasta" | sha256sum -c -
