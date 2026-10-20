#!/usr/bin/env python3
"""Variant table for VariantDatabase from UniProt humsavar, with donor-genotype truth from the PXD022787
personal protein FASTA.

A variant is 'carried' if the variant window (the substitution with up to 6 flanking residues each side,
I = L) occurs in a personal protein sequence that a donor haplotype carries (header tag with xHy or A),
'absent' if only the reference window occurs there, else 'unlabelled'.
Writes variants.tsv (accession, position, ref, alt, id) and variants_tiers.tsv (+ tier: sample if carried,
else curated), and variant_truth.tsv."""
import gzip, re, sys
from collections import defaultdict

AA3 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
       'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}

def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    name, seq = None, []
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq)
                name, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if name is not None: yield name, ''.join(seq)

def carried_sequences(personal):
    """sequences a donor haplotype carries (I->L)"""
    out = []
    for header, seq in read_fasta(personal):
        tag = header.split()[0]
        parts = tag.split('_')[1:]
        if any(p == 'A' or re.fullmatch(r'\dH\d', p) for p in parts):
            out.append(seq.replace('I', 'L'))
    return out

def main(sprot, humsavar, personal, outdir, flank=6):
    prot = {}
    for header, seq in read_fasta(sprot):
        acc = header.split()[0].split('|')[1]
        prot[acc] = seq
    variants = []
    pat = re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')
    with open(humsavar) as f:
        for line in f:
            fields = line.split()
            if len(fields) < 5 or not fields[2].startswith('VAR_'): continue
            m = pat.match(fields[3])
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3: continue
            acc, pos, ref, alt = fields[1], int(m.group(2)), AA3[m.group(1)], AA3[m.group(3)]
            seq = prot.get(acc)
            if seq is None or pos > len(seq) or seq[pos - 1] != ref or ref == alt: continue
            variants.append((acc, pos, ref, alt, fields[2], fields[4]))
    # windows to look up
    need = {}
    for acc, pos, ref, alt, vid, cat in variants:
        seq = prot[acc]
        a, b = max(0, pos - 1 - flank), min(len(seq), pos + flank)
        w_ref = seq[a:b].replace('I', 'L')
        w_alt = (seq[a:pos - 1] + alt + seq[pos:b]).replace('I', 'L')
        need[(acc, pos, alt)] = (w_ref, w_alt)
    lengths = {len(w) for pair in need.values() for w in pair}
    wanted = {w for pair in need.values() for w in pair}
    found = set()
    for s in carried_sequences(personal):
        for L in lengths:
            for i in range(len(s) - L + 1):
                k = s[i:i + L]
                if k in wanted: found.add(k)
    n = defaultdict(int)
    with open(f'{outdir}/variants.tsv', 'w') as v1, open(f'{outdir}/variants_tiers.tsv', 'w') as v2, \
         open(f'{outdir}/variant_truth.tsv', 'w') as vt:
        v1.write('accession\tposition\tref\talt\tid\n')
        v2.write('accession\tposition\tref\talt\tid\ttier\n')
        vt.write('accession\tposition\tref\talt\tid\tcategory\ttruth\n')
        seen = set()
        for acc, pos, ref, alt, vid, cat in variants:
            if (acc, pos, alt) in seen: continue
            seen.add((acc, pos, alt))
            w_ref, w_alt = need[(acc, pos, alt)]
            if w_ref == w_alt: truth = 'IL'  # I<->L: indistinguishable
            elif w_alt in found: truth = 'carried'
            elif w_ref in found: truth = 'absent'
            else: truth = 'unlabelled'
            n[truth] += 1
            v1.write(f'{acc}\t{pos}\t{ref}\t{alt}\t{vid}\n')
            v2.write(f'{acc}\t{pos}\t{ref}\t{alt}\t{vid}\t{"sample" if truth == "carried" else "curated"}\n')
            vt.write(f'{acc}\t{pos}\t{ref}\t{alt}\t{vid}\t{cat}\t{truth}\n')
    print(len(seen), 'variants;', dict(n), file=sys.stderr)

if __name__ == '__main__':
    main(*sys.argv[1:5])
