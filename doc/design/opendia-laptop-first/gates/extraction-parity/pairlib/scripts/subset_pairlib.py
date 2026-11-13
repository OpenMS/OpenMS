"""Pair-preserving 100K library (PLAN_PAIRLIB, label-free construction, run on codon with python3.9 stdlib).

Human side: every non-entrapment entry of the ProteoBench 100K sampled FASTA (unchanged).
Entrapment side: for each of those entries, its own FDRBench shuffle (`_p_target`) from the full ProteoBench FASTA
(entries come in consecutive target / p_target pairs; checked here). The 100K FASTA's own entrapments are dropped.
Library rows (targets and their DECOY_ rows) are copied verbatim from the full predicted library
(`full_fasta/...entrapment.decoys.tsv.zst`, the same prediction pipeline as the 100K library). ProteinId lists are
cut to the selected accessions (shared peptides) whose decoy status matches the row's Decoy flag; rows left without
an accession are dropped. Writes the subset TSV, the paired FASTA and a summary.
Usage: zstd -dc FULL.tsv.zst | python3 subset_pairlib.py FASTA100K FASTAFULL OUT.tsv OUT.fasta > summary.txt"""
import sys

F100, FFULL, OUT, OUTF = sys.argv[1:5]


def read_fasta(path):
    heads, seqs = [], []
    for line in open(path):
        if line[0] == '>':
            heads.append(line[1:].strip())
        else:
            seqs.append(line.strip())
    return heads, seqs


def acc(head):
    return head.split('|')[1]


h100, s100 = read_fasta(F100)
hum = {acc(h): s for h, s in zip(h100, s100) if not acc(h).endswith('_p_target')}
hf, sf = read_fasta(FFULL)
assert len(hf) % 2 == 0
pair, fullhead = {}, {}
for i in range(0, len(hf), 2):
    a, b = acc(hf[i]), acc(hf[i + 1])
    assert a.endswith('_target') and not a.endswith('_p_target') and b.endswith('_p_target'), (a, b)
    assert sorted(sf[i]) == sorted(sf[i + 1]), (a, b)
    pair[a] = (b, sf[i + 1]); fullhead[a] = hf[i]; fullhead[b] = hf[i + 1]
missing = [a for a in hum if a not in pair]
sel = set(hum) | {pair[a][0] for a in hum if a in pair}
with open(OUTF, 'w') as fo:
    for a in sorted(hum):
        if a in pair:
            fo.write('>' + fullhead[a] + '\n' + hum[a] + '\n>' + fullhead[pair[a][0]] + '\n' + pair[a][1] + '\n')

hdr = sys.stdin.readline()
cols = hdr.rstrip('\n').split('\t')
ip = cols.index('ProteinId'); idec = cols.index('Decoy'); ig = cols.index('TransitionGroupId')
n_in = n_out = 0
groups = {0: set(), 1: set()}
seen_acc = set()
with open(OUT, 'w') as fo:
    fo.write(hdr)
    for line in sys.stdin:
        n_in += 1
        f = line.rstrip('\n').split('\t')
        keep = []
        row_decoy = f[idec] == '1'
        for p in f[ip].split(';'):
            # The full library merges protein lists of a peptide that is both a target and another target's decoy;
            # keep only entries whose decoy status matches the row.
            if p.startswith('DECOY_') != row_decoy:
                continue
            a = p[6:] if p.startswith('DECOY_') else p
            parts = a.split('|')
            if len(parts) > 1 and parts[1] in sel:
                keep.append(p); seen_acc.add(parts[1])
        if not keep:
            continue
        f[ip] = ';'.join(keep)
        fo.write('\t'.join(f) + '\n')
        n_out += 1
        groups[int(f[idec])].add(f[ig])
print('100K human entries', len(hum), 'without a pair in the full FASTA', len(missing), missing[:5])
print('selected accessions', len(sel), 'found in library', len(seen_acc))
print('full library rows', n_in, 'kept', n_out)
print('kept target precursors', len(groups[0]), 'decoy precursors', len(groups[1]))
