"""(a) DIA-NN 2.2 on our exact library (raw01): do its accepted entrapments also sit on co-eluting same-composition
reference peptides? (b) For our accepted entrapment winners with such a parent: how many of the 6 library fragments
coincide (10 ppm) with the parent's theoretical b/y ions (charge 1-2)? Evaluation only."""
import sys, numpy as np, pandas as pd
sys.argv = [sys.argv[0], 'D11FWT', '../xicids-20261112', 'xic']
exec(open('diagnose.py').read().split("for raw in ('raw01', 'raw02'):")[0])
MONO = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=160.03065, L=113.08406, I=113.08406,
            N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
            R=156.10111, Y=163.06333, W=186.07931)
P, H2O = 1.007276, 18.010565
def ions(seq):
    m = [MONO[a] for a in seq]; out = []
    for i in range(1, len(seq)):
        b = sum(m[:i]); y = sum(m[i:]) + H2O
        for z in (1, 2):
            out += [(b + z * P) / z, (y + z * P) / z]
    return np.sort(np.array(out))
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0) & ref.Run.str.endswith('_01')].copy()
comp = lambda s: ''.join(sorted(s))
ref['key'] = ref['Stripped.Sequence'].map(comp) + '/' + ref['Precursor.Charge'].astype(str)
par = ref.groupby('key').agg(prt=('RT', lambda x: list(x * 60)), pseq=('Stripped.Sequence', lambda x: sorted(set(x))))

# (a) DIA-NN on our library.
d = pd.read_parquet('../diann-ourlib-20261112/report_raw01.parquet', columns=['Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge', 'RT', 'Q.Value', 'Decoy'])
lab = pd.read_sql("""select pe.MODIFIED_SEQUENCE||p.CHARGE pid, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=m.PROTEIN_ID where p.DECOY=0 group by pid""", c).set_index('pid').e
d = d[d.Decoy == 0].copy(); d['e'] = d['Precursor.Id'].map(lab)
d['key'] = d['Stripped.Sequence'].map(comp) + '/' + d['Precursor.Charge'].astype(str)
d = d.join(par, on='key')
d['other_parent'] = [int(isinstance(r.prt, list) and any(s != r._2 for s in r.pseq) and any(abs(t - r.RT * 60) <= 10 for t in r.prt))
                     for r in d.itertuples()]
print('(a) DIA-NN 2.2 on our library, raw01: accepted targets with a different same-composition reference peptide within 10 s')
for q in (0.01, 0.05, 0.10):
    a = d[d['Q.Value'] <= q]
    print(f'   q <= {q:.2f}: human {a[a.e == 0].other_parent.sum()}/{(a.e == 0).sum()}  entrap {a[a.e == 1].other_parent.sum()}/{(a.e == 1).sum()}'
          f'  -> FDP_c {a.e.sum() * (1 + 1 / R) / len(a):.4f}, without parent-explained entrapments '
          f'{(a.e.sum() - a[a.e == 1].other_parent.sum()) * (1 + 1 / R) / len(a):.4f}')

# (b) Our winners: fragment coincidence with the parent.
seq = pd.read_sql("""select p.ID precursor_id, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID""", c).drop_duplicates('precursor_id').set_index('precursor_id').seq
frag = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz from TRANSITION_PRECURSOR_MAPPING m join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
pop = pd.read_csv(f'{ROOT}/raw01/population.tsv', sep='\t')
s = np.load(f'{ROOT}/raw01/scores_{ARM}.npy').astype(np.float64)
w, best = tdc(pop, s)
best = best[best.is_decoy == 0].join(pre[['charge']], on='precursor_id')
best['e'] = best.precursor_id.map(ent); best['seq'] = best.precursor_id.map(seq)
best['key'] = best.seq.map(comp) + '/' + best.charge.astype(str)
cut10 = best[best.q <= 0.10].s.min()
top = best[best.s >= cut10].join(par, on='key')
fr = frag[frag.precursor_id.isin(top.precursor_id)].groupby('precursor_id').mz.apply(np.array)
print('(b) our raw01 winners above the q <= 0.10 cut with a co-eluting (10 s) different same-composition reference peptide:')
for r in top.itertuples():
    if not isinstance(r.prt, list):
        continue
    others = [p for p, t in zip(r.pseq, [r.prt] * len(r.pseq)) if p != r.seq]
    if not others or not any(abs(t - r.rt_apex) <= 10 for t in r.prt):
        continue
    th = ions(others[0]); lib = fr[r.precursor_id]
    shared = int(sum(np.min(np.abs(th - m)) <= m * 10e-6 for m in lib))
    print(f'   {"ENTRAP" if r.e == 1 else "human "} {r.seq}/{r.charge} parent {others[0]} q {r.q:.4f}: {shared}/{len(lib)} library fragments are parent ions')
