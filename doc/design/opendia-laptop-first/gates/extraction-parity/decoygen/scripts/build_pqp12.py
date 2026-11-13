"""PLAN_DECOYGEN data step (replaces PAIRLIB's build_pqp12.py in the arm runs; run under the arm's bind at the PAIRLIB path):
12-fragment PQP for the voxel populations, used only for the -out_xipm export.
Precursors (IDs, m/z, RT, decoy flags, peptides) are copied unchanged from the arm's library. Transitions are replaced:
  - a target gets its DIALibGen 0.13.1 top-12 b/y fragments (lib/pred12.parquet, the PAIRLIB prediction, unchanged);
  - a decoy gets its OWN DIALibGen 0.13.1 top-12 prediction (lib/decoy12.tsv, written by build_lib.py from decoys_p12.parquet,
    the same settings as pred12). This is the one change from PAIRLIB, where a decoy copied its target's annotations and
    intensities. A decoy without a 12-fragment prediction keeps its 6 library transitions (count printed).
Fragment m/z is recomputed on each precursor's own modified sequence; checks: the mass code reproduces every library
fragment m/z within 2e-3, and the predicted m/z of every target and decoy fragment within 2e-3. Only population precursors
are kept."""
import re, shutil, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
V = B + '/pairlib-20261113'  # bound to decoygen-20261113/arm<X> on codon
SRC = V + '/lib/library.pqp'
OUT = V + '/lib12/library12_pop.pqp'
PROTON, H2O = 1.007276467, 18.010564684
AA = dict(G=57.021464, A=71.037114, S=87.032028, P=97.052764, V=99.068414, T=101.047679, C=103.009185, L=113.084064,
          I=113.084064, N=114.042927, D=115.026943, Q=128.058578, K=128.094963, E=129.042593, M=131.040485, H=137.058912,
          F=147.068414, R=156.101111, Y=163.06332, W=186.079313)
MOD = {'UniMod:4': 57.021464, 'UniMod:35': 15.994915}


def residues(modseq):
    return np.array([AA[a] + (MOD[m] if m else 0.0) for a, m in re.findall(r'([A-Z])(?:\((UniMod:\d+)\))?', modseq)])


def frag_mz(res, ion, ordinal, z):
    m = res[:ordinal].sum() if ion == 'b' else res[-ordinal:].sum() + H2O
    return (m + z * PROTON) / z


def norm(s):
    return s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')


src = sqlite3.connect('file:' + SRC + '?mode=ro', uri=True)
prec = pd.read_sql("""select p.ID id, p.TRAML_ID traml, p.CHARGE z, p.DECOY decoy, pe.MODIFIED_SEQUENCE ms from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", src)
tr = pd.read_sql("""select m.PRECURSOR_ID id, t.PRODUCT_MZ mz, t.TYPE ion, t.ORDINAL ord, t.CHARGE fz, t.LIBRARY_INTENSITY li
  from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID""", src)
seq_of = dict(zip(prec.id, prec.ms)); cache = {}
err = np.empty(len(tr))
for k, (i, mz, ion, o, fz) in enumerate(tr[['id', 'mz', 'ion', 'ord', 'fz']].itertuples(index=False)):
    r = cache.get(i)
    if r is None: r = cache[i] = residues(seq_of[i])
    err[k] = abs(frag_mz(r, ion, o, fz) - mz)
print('library fragment m/z check: max abs error', err.max(), 'n', len(err), flush=True)
assert err.max() < 2e-3

need = set()
for raw in ('raw01', 'raw02'):
    need |= set(pd.read_csv(f'{V}/{raw}/population.tsv', sep='\t', usecols=['precursor_id']).precursor_id)
prec = prec[prec.id.isin(need)].copy()

pred = pq.read_table(V + '/lib/pred12.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Product.Mz', 'Relative.Intensity',
                                                        'Fragment.Type', 'Fragment.Charge', 'Fragment.Series.Number']).to_pandas()
pred['key'] = pred['Modified.Sequence'].map(norm) + '/' + pred['Precursor.Charge'].astype(str)
pred = pred.drop_duplicates('key').set_index('key')
d12 = pd.read_csv(V + '/lib/decoy12.tsv', sep='\t')
d12 = {t: g.sort_values('rank') for t, g in d12.groupby('traml')}
lib6 = {i: g for i, g in tr[tr.id.isin(prec.id)].groupby('id')}

rows = []; tid = 0; dev = []; fallback = 0; miss_t = 0
for p in prec.itertuples(index=False):
    res = residues(p.ms)
    if p.decoy == 0:
        key = p.ms + '/' + str(p.z)
        if key not in pred.index:
            miss_t += 1; continue
        q = pred.loc[key]
        frags = [({2: 'b', 5: 'y'}[int(ft)], int(o), int(fz), float(mz), float(li)) for mz, li, ft, fz, o in
                 zip(q['Product.Mz'], q['Relative.Intensity'], q['Fragment.Type'], q['Fragment.Charge'], q['Fragment.Series.Number'])]
    elif p.traml in d12:
        g = d12[p.traml]
        frags = list(zip(g.ion, g.ord.astype(int), g.fz.astype(int), g.mz.astype(float), g.li.astype(float)))
    else:
        fallback += 1; g = lib6[p.id]
        frags = list(zip(g.ion, g.ord.astype(int), g.fz.astype(int), g.mz.astype(float), g.li.astype(float)))
    for ion, o, fz, mz, li in frags:
        m = frag_mz(res, ion, o, fz); dev.append(abs(m - mz))
        rows.append((tid, f'{p.traml}_{ion}{o}^{fz}', m, fz, ion, f'{ion}{o}^{fz}', o, 1, 0, 1, li, int(p.decoy), int(p.id))); tid += 1
print('targets without pred12 (no transitions)', miss_t, '| decoys without a 12-fragment prediction (6 library transitions kept)', fallback)
print('predicted fragment m/z against own sequence: max abs error', max(dev), flush=True)
assert max(dev) < 2e-3
t = pd.DataFrame(rows, columns=['ID', 'TRAML_ID', 'PRODUCT_MZ', 'CHARGE', 'TYPE', 'ANNOTATION', 'ORDINAL', 'DETECTING',
                                'IDENTIFYING', 'QUANTIFYING', 'LIBRARY_INTENSITY', 'DECOY', 'PRECURSOR_ID'])
shutil.copyfile(SRC, OUT)
c = sqlite3.connect(OUT)
pep = pd.read_sql('select PRECURSOR_ID, PEPTIDE_ID from PRECURSOR_PEPTIDE_MAPPING', c).set_index('PRECURSOR_ID').PEPTIDE_ID
keep = ','.join(map(str, sorted(prec.id)))
for tab in ('TRANSITION', 'TRANSITION_PRECURSOR_MAPPING', 'TRANSITION_PEPTIDE_MAPPING'):
    c.execute(f'delete from {tab}')
c.execute(f'delete from PRECURSOR where ID not in ({keep})')
c.execute('delete from PRECURSOR_PEPTIDE_MAPPING where PRECURSOR_ID not in (select ID from PRECURSOR)')
c.executemany('insert into TRANSITION values (?,?,?,?,?,?,?,?,?,?,?,?)', t.iloc[:, :12].itertuples(index=False, name=None))
c.executemany('insert into TRANSITION_PRECURSOR_MAPPING values (?,?)', zip(t.ID.tolist(), t.PRECURSOR_ID.tolist()))
c.executemany('insert into TRANSITION_PEPTIDE_MAPPING values (?,?)', zip(t.ID.tolist(), pep.loc[t.PRECURSOR_ID].tolist()))
c.commit(); c.execute('vacuum'); c.close()
t[['PRECURSOR_ID', 'ANNOTATION', 'PRODUCT_MZ', 'LIBRARY_INTENSITY', 'DECOY']].to_csv(OUT.replace('.pqp', '.transitions.tsv'), sep='\t', index=False)
print('precursors', len(prec), 'targets', int((prec.decoy == 0).sum()), 'decoys', int((prec.decoy == 1).sum()), 'transitions', len(t),
      'per precursor', t.groupby('PRECURSOR_ID').size().value_counts().sort_index().to_dict())
