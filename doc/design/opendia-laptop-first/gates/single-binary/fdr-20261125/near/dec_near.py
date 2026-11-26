#!/usr/bin/env python3
"""Top-region decoys vs entrapments by near-human class. Step 1 (dump): write unmodified sequences of decoys and of
entrapment targets with tool q <= QMAX in the G0 arm of a set. Step 2 (report): given hd tables from strict_list.py
SEQS mode, print the class shares for decoys and entrapments, by q band, against the library controls.
Usage: dec_near.py dump <lib.pqp> <runs dir> <SET> <arm> <out prefix>   env RELABEL
       dec_near.py report <out prefix> <hd all.tsv> <lib decoy near share> <lib entrapment near share>"""
import sqlite3, sys, glob, os, pandas as pd
QB = [(0, 0.01), (0.01, 0.05), (0.05, 0.3)]
if sys.argv[1] == 'dump':
    L, R, S, arm, out = sys.argv[2:7]
    c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
    pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
      min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
      from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
      join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
    if os.environ.get('RELABEL'):
        rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); pq.loc[(pq.d == 0) & pq.seq.isin(rl), 'e'] = 0
    pi = pq.set_index('precursor_id'); rows = []
    for d in sorted(glob.glob(f'{R}/{S}_{arm}_s*')):
        f = [x for x in glob.glob(f'{d}/*.precursors.tsv*') if '/precursors.tsv' not in x][0]
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[(t.score > -1e29) & (t.q_value <= QB[-1][1])].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        t = t.join(pi[['seq', 'e']], on='precursor_id')
        t['cls'] = (t.is_decoy == 1).map({True: 'D', False: None}); t.loc[(t.is_decoy == 0) & (t.e == 1), 'cls'] = 'E'
        rows.append(t[t.cls.notna()][['seq', 'cls', 'q_value']].assign(seed=d[-1]))
    T = pd.concat(rows); T.to_csv(out + '.rows.tsv', sep='\t', index=False)
    open(out + '.seqs', 'w').write('\n'.join(sorted(set(T.seq))) + '\n')
    print(f'{S} {arm}: rows {len(T)} (decoy {int((T.cls == "D").sum())}, entrapment {int((T.cls == "E").sum())}), unique seqs {T.seq.nunique()}')
else:
    out, hdf, ld, le = sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5])
    T = pd.read_csv(out + '.rows.tsv', sep='\t').join(pd.read_csv(hdf, sep='\t').set_index('seq').hd, on='seq')
    print(f'library near (hd <= 2) share: decoys {ld:.3f}, entrapments {le:.3f}')
    for a, b in QB:
        g = T[(T.q_value > a) & (T.q_value <= b)]
        s = ' | '.join(f'{k} n {len(h)} hd1 {(h.hd == 1).mean():.3f} hd<=2 {(h.hd <= 2).mean():.3f} (x{(h.hd <= 2).mean() / (ld if k == "D" else le):.2f})'
                       for k, h in g.groupby('cls'))
        print(f'  q ({a}, {b}]: {s}')
