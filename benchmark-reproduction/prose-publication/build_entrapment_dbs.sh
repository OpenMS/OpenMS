#!/usr/bin/env bash
# Builds the target/entrapment/decoy databases of the ProSE publication benchmark (DATASETS.md, section 3) with
# DecoyDatabase from OpenMS develop (>= 916a7d5, PR #10629).
#
# usage: build_entrapment_dbs.sh <OpenMS bin dir> <work dir>
#
# <work dir>/data must hold the target-only FASTAs hy_targets.fasta, hye_targets.fasta and human_targets.fasta
# (prepare_inputs.sh writes them) and arath_sp.fasta. Outputs go to <work dir>/db, with SHA256SUMS.
set -euo pipefail
BIN=$1
W=$2
D=$W/data
O=$W/db
mkdir -p "$O"

# Decoys: full protein reversal with DECOY_ prefix, as the #10364 suite (keepPeptideNTerm/CTerm false turns off the
# peptide-wise reversal). Entrapment: peptide-wise shuffle, C-terminal residue and K, R and P fixed (DecoyDatabase's only
# design). The enzyme is the one of the search, Trypsin/P.
decoy() {
  "$BIN/DecoyDatabase" -in "$1" -out "$2" -enzyme Trypsin/P \
    -Decoy:keepPeptideNTerm false -Decoy:keepPeptideCTerm false "${@:3}" -no_progress
}

# 1. Shuffled entrapment, r = 1, three realisations per database
for db in hy hye human; do
  for seed in 1 2 3; do
    decoy "$D/${db}_targets.fasta" "$O/${db}_entrap_s${seed}_td.fasta" -entrapment_ratio 1 -seed "$seed"
  done
done

# 2. Velos natural entrapment: the sample is yeast + the 48 UPS1 proteins, so every human protein other than UPS1 (and
# not a contaminant) gets the ENTRAP_ prefix (as make_velos_natural.py of #10629 does), then decoys as above
python3 -I - "$D/hy_targets.fasta" "$O/hy_velos_natural_targets.fasta" <<'EOF'
import sys
UPS1 = set("""P00915 P00918 P01031 P69905 P68871 P41159 P02768 P0CG48 P0CG47 P62987 P62979 P04040 P00167 P01133 P02144 P15559
P62937 Q06830 P63165 P00709 P06730 P12081 P61626 Q15843 P02753 P16083 P63279 P01008 P61769 P55957 O76070 P08263 P01344
P01127 P10599 P99999 P06396 P09211 P01112 P01579 P02787 P02788 P05413 P10636 P10145 P02741 O00762 P01375 P51965 P08758""".split())
n = 0
with open(sys.argv[2], "w") as out:
    for line in open(sys.argv[1]):
        if line.startswith(">"):
            h = line[1:]
            acc = h.split("|")[1] if "|" in h else ""
            if "_HUMAN" in h.split()[0] and "Cont_" not in h and acc not in UPS1:
                line = ">ENTRAP_" + h
                n += 1
        out.write(line)
print(sys.argv[2], "entrapment proteins:", n)
EOF
decoy "$O/hy_velos_natural_targets.fasta" "$O/hy_velos_natural_td.fasta"
rm "$O/hy_velos_natural_targets.fasta"

# 3. Foreign-proteome entrapment (reference for the shuffle): reviewed Arabidopsis thaliana (UP000006548) as ENTRAP_
# targets next to the HYE and human targets
for db in hye human; do
  { cat "$D/${db}_targets.fasta"; sed 's/^>/>ENTRAP_/' "$D/arath_sp.fasta"; } > "$O/${db}_arath_targets.fasta"
  decoy "$O/${db}_arath_targets.fasta" "$O/${db}_arath_td.fasta"
  rm "$O/${db}_arath_targets.fasta"
done

(cd "$O" && sha256sum ./*_td.fasta > SHA256SUMS)
cat "$O/SHA256SUMS"
