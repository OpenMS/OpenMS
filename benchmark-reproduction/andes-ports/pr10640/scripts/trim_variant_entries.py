#!/usr/bin/env python3
"""Shorten the variant and knockoff entries of a VariantDatabase output for searching.

VariantDatabase writes every variant and knockoff as a full-length protein copy (247 M residues for 83,473
UniProt variants + 65,044 knockoffs, 22x Swiss-Prot human; TTN alone 47.6 M), which does not fit a 15 GB search.
Each entry is cut to a window around its substitution that starts and ends at Trypsin/P cleavage sites, at least
4 cleavage sites and 45 residues from the site on each side (and not starting with M), so that every peptide
with <= 2 missed cleavages that holds the site or ends/starts next to it is kept; the other peptides of an entry
are canonical ones, which the canonical protein provides anyway. Accessions and descriptions are unchanged.
Usage: trim_variant_entries.py in.fasta out.fasta"""
import re, sys

VAR = re.compile(r'_VAR_[A-Z](\d+)[A-Z]$')
KO = re.compile(r'_KO\d+_[A-Z](\d+)[A-Z]_VAR_[A-Z]\d+[A-Z]$')

def window(seq, p, sites=4, min_res=45):
    cuts = [i + 1 for i, c in enumerate(seq[:-1]) if c in 'KR']  # boundaries (start index of next segment)
    left = [c for c in cuts if c <= p]
    a = left[-sites] if len(left) >= sites else 0
    a = min(a, max(0, p - min_res))
    a = max([c for c in [0] + cuts if c <= a])  # snap to a boundary
    while a > 0 and seq[a] == 'M':
        a = max([c for c in [0] + cuts if c < a])
    right = [c for c in cuts if c > p + 1]
    b = right[sites - 1] if len(right) >= sites else len(seq)
    b = max(b, min(len(seq), p + 1 + min_res))
    later = [c for c in cuts if c >= b]
    b = later[0] if later else len(seq)
    return a, b

def main(inp, out):
    n = before = after = 0
    with open(inp) as f, open(out, 'w') as o:
        header, seq = None, []
        def flush():
            nonlocal n, before, after
            if header is None: return
            s = ''.join(seq)
            acc = header[1:].split()[0]
            m = KO.search(acc) or VAR.search(acc)
            before += len(s)
            if m:
                a, b = window(s, int(m.group(1)) - 1)
                s = s[a:b]; n += 1
            after += len(s)
            o.write(header + '\n')
            for i in range(0, len(s), 60): o.write(s[i:i + 60] + '\n')
        for line in f:
            if line.startswith('>'):
                flush(); header, seq = line.rstrip('\n'), []
            else:
                seq.append(line.strip())
        flush()
    print(f'{n} entries trimmed; residues {before/1e6:.1f} M -> {after/1e6:.1f} M', file=sys.stderr)

if __name__ == '__main__':
    main(*sys.argv[1:3])
