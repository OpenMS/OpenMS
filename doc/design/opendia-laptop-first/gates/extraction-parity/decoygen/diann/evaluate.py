"""DIA-NN evaluation on the PAIRLIB library family (evaluation only; reads entrapment labels from the library).
Per report: precursor-level target IDs at DIA-NN's run q (Q.Value) <= 0.01 and 0.05, entrapment FDP_c = E (1 + 1/r) / N
(r = entrapment / human target precursors of the library, as evaluate_v2.py), N at FDP_c 1 % and 2 % (targets ranked
by Q.Value then CScore), and whether library decoys were used (decoy rows in the report keyed to library decoy
precursors). Entrapment label: the precursor's library protein accessions all carry _p_target (DIA-NN reports drop it,
so labels are mapped through the library by modified sequence + charge, as in diann-pairlib-20261113).
Usage: evaluate.py <tag>/<raw> ..."""
import sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
lib = pd.read_csv(B + '/pairlib-20261113/lib/pairlib.decoys.tsv', sep='\t',
                  usecols=['ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId', 'Decoy'], dtype=str)
lib = lib[lib.Decoy == '0'].drop_duplicates(['ModifiedPeptideSequence', 'PrecursorCharge'])
lib['ent'] = lib.ProteinId.str.split(';').map(lambda ps: all('_p_target' in p for p in ps))
lib['key'] = lib.ModifiedPeptideSequence + '/' + lib.PrecursorCharge
ent = lib.set_index('key').ent
r = ent.sum() / (~ent).sum()
k = 1 + 1 / r
print(f'library targets {len(ent)}, entrapment {int(ent.sum())}, r {r:.4f}')


def n_at(sorted_ent, t):
    f = np.cumsum(sorted_ent) * k / np.arange(1, len(sorted_ent) + 1)
    ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0


for path in sys.argv[1:]:
    rep = pd.read_parquet(f'{B}/diann-decoys-20261114/{path}/report.parquet' if not path.startswith('/') else path)
    dec = rep['Decoy'].sum() if 'Decoy' in rep else np.nan
    t = rep[rep.get('Decoy', 0) == 0].copy()
    t['key'] = t['Modified.Sequence'] + '/' + t['Precursor.Charge'].astype(str)
    t = t.sort_values(['Q.Value', 'CScore'] if 'CScore' in t else ['Q.Value'], ascending=[True, False] if 'CScore' in t else True)
    t = t.drop_duplicates('key')
    t['ent'] = t.key.map(ent)
    unmapped = t.ent.isna().sum()
    t = t.dropna(subset=['ent'])
    e = t.ent.values.astype(bool)
    line = [f'{path:18s}']
    for q in (0.01, 0.05):
        s = t[t['Q.Value'] <= q]
        line.append(f'q<={q}: N {len(s)} E {int(s.ent.sum())} FDP_c {s.ent.sum() * k / max(len(s), 1):.4f}')
    line.append(f'N@FDP_c 1%/2%: {n_at(e, 0.01)}/{n_at(e, 0.02)}')
    line.append(f'decoy rows {int(dec) if dec == dec else "n/a"}, unmapped targets {unmapped}')
    print(' | '.join(line))
