"""Wrong-peak rate of PIP-ECHO transfers from pseudo-transfers (PIPECHO_PSEUDO dump): for each hidden identification
(precursor K hidden in run m), an accepted transfer of K into m is right if its acceptor is the hidden feature
(acceptor_pseudo_key == K) and a wrong peak otherwise. Reports, per transfer q threshold, the accepted pseudo-transfers,
the observed wrong-peak share, and the recovery (right transfers / hidden identifications).
With the feature removed (PIPECHO_PSEUDO_MODE=absent) every accepted transfer into a hidden pair is wrong; then
"wrong_peak_pct" is 100 and "false_acceptance_pct" (wrong transfers / hidden identifications) is the measure.
Usage: pseudo_eval.py hidden.tsv [label]"""
import sys, pandas as pd
h = pd.read_csv(sys.argv[1], sep='\t'); a = pd.read_csv(sys.argv[1] + '.accepted.tsv', sep='\t')
label = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
hidden = set(zip(h.key, h.hidden_map))
a = a[[(k, m) in hidden for k, m in zip(a.donor_key, a.acceptor_map)]].copy()
a['right'] = a.acceptor_pseudo_key == a.donor_key
print(f'{label}: hidden identifications {len(h)}')
print('q\taccepted\twrong_peak\twrong_peak_pct\trecovery_pct\tfalse_acceptance_pct')
for q in (0.005, 0.01, 0.02, 0.03, 0.05):
    x = a[a.q <= q]
    w = int((~x.right).sum())
    print(f'{q}\t{len(x)}\t{w}\t{100 * w / max(len(x), 1):.2f}\t{100 * x.right.sum() / len(h):.1f}\t{100 * w / len(h):.2f}')
