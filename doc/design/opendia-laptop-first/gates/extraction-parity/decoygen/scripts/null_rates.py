"""PLAN_DECOYGEN label-free pre-check (recorded; no labels): random-match rates of target and decoy windows in null windows.
Null window: a population window whose |apex - calibrated library RT| > 3 sigma (no peptide expected there).
Per class: mean transition_support, share with support >= 4, mean f12_libcos, share f12_libcos > 0.8.
Label-free: target vs decoy only. With --labels the targets are split into human (H) and entrapment (E) and the decoys by
the class of their target (D(H), D(E)); that part is evaluation only.
Usage (run with the arm's root bound to the PAIRLIB path): null_rates.py [--labels]"""
import json, sqlite3, sys, numpy as np, pandas as pd
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113'
LABELS = '--labels' in sys.argv
c = sqlite3.connect(f'file:{O}/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
if LABELS:
    ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID
      group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
    tid = pre[pre.decoy == 0].set_index('traml').precursor_id
    pre['partner'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tid), pre.precursor_id)
    pre['cls'] = np.where(pre.decoy == 1, np.where(pre.partner.map(ent) == 1, 'D(E)', 'D(H)'),
                          np.where(pre.precursor_id.map(ent) == 1, 'E', 'H'))
else:
    pre['cls'] = np.where(pre.decoy == 1, 'D', 'T')
cls = pre.set_index('precursor_id').cls
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t', usecols=['precursor_id'])
    H = np.load(f'{O}/{raw}/heuristics.npy', mmap_mode='r'); F = np.load(f'{O}/{raw}/feat12.npy', mmap_mode='r')
    sig = json.load(open(f'{O}/{raw}/cal/cal.json'))['rt']['sigma_s']
    null = np.asarray(H[:, 3]) > 3 * sig
    sup = np.asarray(H[:, 0])[null]; lc = np.asarray(F[:, 2])[null]; k = pop.precursor_id.values[null]
    d = pd.DataFrame({'cls': cls.loc[k].values, 'sup': sup, 'lc': lc})
    print(f'==== {raw}: null windows {null.sum()} (|RT err| > 3 sigma = {3 * sig:.1f} s)')
    print('  class        n   mean support  support>=4  mean libcos  libcos>0.8')
    for name, g in d.groupby('cls'):
        print(f'  {name:5s} {len(g):9d}   {g.sup.mean():10.3f}  {np.mean(g.sup >= 4):10.4f}  {g.lc.mean():11.3f}  {np.mean(g.lc > 0.8):10.4f}')
    if not LABELS:
        t, dd = d[d.cls == 'T'], d[d.cls == 'D']
        print(f'  decoy / target: support>=4 {np.mean(dd.sup >= 4) / np.mean(t.sup >= 4):.3f}  libcos>0.8 '
              f'{np.mean(dd.lc > 0.8) / np.mean(t.lc > 0.8):.3f}  mean support {dd.sup.mean() / t.sup.mean():.3f}')
