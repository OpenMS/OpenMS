"""DIA-NN 2.2 on the r~7 entrapment library, scored with the same labels, r and bounds as eval_pdec.py.
Usage: eval_diann_r5.py <library.pqp> <report.parquet> [<opendia precursors.tsv.gz> ...]
Entrapment labels are read here only. DIA-NN ran with --qvalue 1, so the report holds every scored target precursor.
At DIA-NN's own Q.Value <= 0.01: N, entrapments E, lower bound E/N, combined bound E(1+1/r)/N, 90 % Poisson interval.
N@FDP_c 1 %: targets ranked by Q.Value then PEP, largest prefix with combined FDP <= 0.01.
Optional OpenDIA cut files (eval_pdec layout) give the overlap of the accepted non-entrapment sets."""
import sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import chi2
L, REP = sys.argv[1:3]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
tg = pd.read_sql("""select p.ID precursor_id, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  where p.DECOY=0 group by p.ID""", c)
tg['key'] = tg.ms + '/' + tg.z.astype(str)
r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
lab = tg.drop_duplicates('key').set_index('key')
d = pd.read_parquet(REP, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'PEP'])
d['key'] = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d.sort_values(['Q.Value', 'PEP']).drop_duplicates('key')
miss = ~d.key.isin(lab.index)
print(f'library: targets {len(tg)}, r = {r:.4f}; DIA-NN report precursors {len(d)}, not matched to the library {int(miss.sum())}')
d = d[~miss].copy(); d['e'] = lab.e.reindex(d.key).values; d['pid'] = lab.precursor_id.reindex(d.key).values
acc = d[d['Q.Value'] <= 0.01]; N, E = len(acc), int(acc.e.sum())
lo, hi = (chi2.ppf(0.05, 2 * E) / 2 if E else 0.0), chi2.ppf(0.95, 2 * E + 2) / 2
fdp = np.cumsum(d.e.values) * k / np.arange(1, len(d) + 1); ok = np.where(fdp <= 0.01)[0]
n1 = int(ok.max() + 1) if len(ok) else 0
print(f'DIA-NN: q<=0.01 N {N}, E {E}, lower {E / N:.4f} (90 % {lo / N:.4f}-{hi / N:.4f}), combined {E * k / N:.4f} '
      f'(90 % {lo * k / N:.4f}-{hi * k / N:.4f}); N@FDP_c 1 % {n1}')
for q in (0.005, 0.02):
    a = d[d['Q.Value'] <= q]; print(f'  at q<={q}: N {len(a)}, E {int(a.e.sum())}, lower {a.e.sum() / len(a):.4f}, combined {a.e.sum() * k / len(a):.4f}')
dn = set(acc[acc.e == 0].pid)
for f in sys.argv[3:]:
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'q_value'])
    t = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]; t = t[t.precursor_id.map(lab.set_index('precursor_id').e) == 0]
    od = set(t.precursor_id)
    print(f'  overlap with {f.split("/")[-1]}: OpenDIA non-entrapment q<=0.01 {len(od)}, in DIA-NN {len(od & dn) / max(len(od), 1):.3f}; '
          f'DIA-NN {len(dn)} in OpenDIA {len(od & dn) / max(len(dn), 1):.3f}')
