"""PLAN_DECOYGEN hypothesis reading: pooled rho over raw01 + raw02 per q cut, from an arm's DIAG_<arm>.txt (diagnose.py output,
section 2 lines 'q <= X cut ...: N n E e D d ...'). Interval as diagnose.py (binomial split of E vs D, 90 %).
r (entrapment / human target precursors) is read from the library as diagnose.py does.
Usage: pooled_rho.py <root> <arm>"""
import re, sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import beta as betad
ROOT, ARM = sys.argv[1:3]
c = sqlite3.connect(f'file:{ROOT}/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
  group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
tgt = pd.read_sql('select ID precursor_id from PRECURSOR where DECOY = 0', c).precursor_id
R = (ent.loc[tgt] == 1).sum() / (ent.loc[tgt] == 0).sum()
k = 1 + 1 / R
tot = {}
txt = open(f'{ROOT}/score/DIAG_{ARM}.txt').read()
for raw_block in txt.split('==== ')[1:]:
    for q, e, d in re.findall(r'q <= ([0-9.]+) cut\s+[-0-9.]+: N\s+\d+ E\s+(\d+) D\s+(\d+)', raw_block.split('3. ')[0]):
        a = tot.setdefault(q, [0, 0, 0]); a[0] += int(e); a[1] += int(d); a[2] += 1
print(f'{ARM} pooled over raws (r = {R:.4f})')
for q, (E, D, n) in tot.items():
    lo = betad.ppf(0.05, E, D + 1) if E > 0 else 0.0
    hi = betad.ppf(0.95, E + 1, D) if D > 0 else 1.0
    f = lambda p: k * p / (1 - p) if p < 1 else np.inf
    print(f'  q <= {q}: raws {n} E {E:4d} D {D:4d} rho {k * E / D if D else np.inf:5.2f} [{f(lo):4.2f}, {f(hi):5.2f}]')
