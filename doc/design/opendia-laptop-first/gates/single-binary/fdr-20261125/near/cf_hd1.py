#!/usr/bin/env python3
"""Decoy-matched Hamming-1 counterfactual (diagnosis). At tool q <= 0.01 per set (arm as in dec_near dump):
the excess of Hamming-1 (one substitution from a human window, I=L) entrapments over what the decoy Hamming-1 share
predicts. E_cf = E_not_hd1 / (1 - share_hd1(accepted decoys)); combined_cf = E_cf (1 + 1/r) / N, with r and N as in
eval_strict.py relabel. Usage: cf_hd1.py <out prefix> <hd all.tsv> <N total> <r>"""
import sys, pandas as pd
out, hdf, N, r = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
T = pd.read_csv(out + '.rows.tsv', sep='\t').join(pd.read_csv(hdf, sep='\t').set_index('seq').hd, on='seq')
a = T[T.q_value <= 0.01]; E, D = a[a.cls == 'E'], a[a.cls == 'D']
se, sd = (E.hd == 1).mean(), (D.hd == 1).mean(); Enon = int((E.hd != 1).sum()); Ecf = Enon / (1 - sd); k = 1 + 1 / r
print(f'{out}: N {N} E {len(E)} (hd1 {se:.3f}) decoys {len(D)} (hd1 {sd:.3f}) | combined {len(E) * k / N:.4f} -> decoy-matched hd1 {Ecf * k / N:.4f}')
