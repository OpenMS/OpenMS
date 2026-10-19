"""Velos natural entrapment (frozen pr-prototypes-20261001/entrapment_velos.py): the sample is yeast + the 48 UPS1
proteins, so every reviewed human protein other than UPS1 (and not a contaminant) is entrapment. EntrapmentEvaluator
recognises entrapment by an accession prefix, so this writes hy_td.fasta with 'ENTRAP_' put in front of those target
accessions. Sequences, order and decoys (DECOY_<original accession>) are unchanged.
usage: make_velos_natural.py hy_td.fasta out.fasta"""
import hashlib, sys
UPS1 = """P00915 P00918 P01031 P69905 P68871 P41159 P02768 P0CG48 P0CG47 P62987 P62979 P04040 P00167 P01133 P02144 P15559
P62937 Q06830 P63165 P00709 P06730 P12081 P61626 Q15843 P02753 P16083 P63279 P01008 P61769 P55957 O76070 P08263 P01344
P01127 P10599 P99999 P06396 P09211 P01112 P01579 P02787 P02788 P05413 P10636 P10145 P02741 O00762 P01375 P51965 P08758""".split()
src, dst = sys.argv[1], sys.argv[2]
n = 0
with open(dst, 'w') as out:
    for line in open(src):
        if line.startswith('>') and not line.startswith('>DECOY_'):
            h = line[1:]
            acc = h.split('|')[1] if '|' in h else ''
            if '_HUMAN' in h.split()[0] and 'Cont_' not in h and acc not in UPS1:
                line = '>ENTRAP_' + h; n += 1
        out.write(line)
print(dst, 'entrapment proteins', n, 'sha256', hashlib.sha256(open(dst, 'rb').read()).hexdigest())
