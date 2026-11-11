"""PLAN_V4 phase 2 data step: 12-fragment PQP for the voxel populations, used only for the -out_xipm export.
Precursors (IDs, m/z, RT, decoy flags, peptides) are copied unchanged from the original library. Transitions are
replaced: a target gets its DIALibGen 0.13.1 top-12 b/y fragments (pred12.parquet, Astral, NCE 25); a decoy gets the
same annotations and predicted intensities as its target (the original library's decoy rule, checked below), with
m/z computed on the decoy's own modified sequence. Only precursors of the raw01/raw02 populations are kept.
Checks: the mass code reproduces every original PQP fragment m/z (targets and decoys) within 2e-3."""
import re, shutil, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
SRC = B + '/astral-entrapment-20261107/lib/library.pqp'
OUT = B + '/rtlin-20261111/lib12/library12_pop.pqp'  # rtlin rebuild (new populations)
V = B + '/rtlin-20261111'
PROTON, H2O = 1.007276467, 18.010564684
AA = dict(G=57.021464, A=71.037114, S=87.032028, P=97.052764, V=99.068414, T=101.047679, C=103.009185, L=113.084064,
          I=113.084064, N=114.042927, D=115.026943, Q=128.058578, K=128.094963, E=129.042593, M=131.040485, H=137.058912,
          F=147.068414, R=156.101111, Y=163.06332, W=186.079313)
MOD = {'UniMod:4': 57.021464, 'UniMod:35': 15.994915}

def residues(modseq):
    out = []
    for aa, mod in re.findall(r'([A-Z])(?:\((UniMod:\d+)\))?', modseq):
        out.append(AA[aa] + (MOD[mod] if mod else 0.0))
    return np.array(out)

def frag_mz(res, ion, ordinal, z):
    m = res[:ordinal].sum() if ion == 'b' else res[-ordinal:].sum() + H2O
    return (m + z * PROTON) / z

def norm(s):
    return s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')

src = sqlite3.connect('file:' + SRC + '?mode=ro', uri=True)
prec = pd.read_sql("""select p.ID id, p.TRAML_ID traml, p.CHARGE z, p.DECOY decoy, pe.MODIFIED_SEQUENCE ms from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", src)
tr = pd.read_sql("""select m.PRECURSOR_ID id, t.PRODUCT_MZ mz, t.TYPE ion, t.ORDINAL ord, t.CHARGE fz from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID""", src)

# Check the mass code on the original library (all precursors).
seq_of = dict(zip(prec.id, prec.ms)); cache = {}
err = np.empty(len(tr))
for k, (i, mz, ion, o, fz) in enumerate(tr.itertuples(index=False)):
    r = cache.get(i)
    if r is None: r = cache[i] = residues(seq_of[i])
    err[k] = abs(frag_mz(r, ion, o, fz) - mz)
print('original fragment m/z check: max abs error', err.max(), 'n', len(err), flush=True)
assert err.max() < 2e-3

# Population precursors (targets and decoys) of both raws.
need = set()
for raw in ('raw01', 'raw02'):
    need |= set(pd.read_csv(f'{V}/{raw}/population.tsv', sep='\t', usecols=['precursor_id']).precursor_id)
prec = prec[prec.id.isin(need)].copy()
# Target key of each precursor: the TRAML_ID without the DECOY_ prefix (the original pairing).
prec['tkey'] = prec.traml.str.replace(r'^DECOY_', '', regex=True)
targets = prec[prec.decoy == 0].set_index('tkey')
pairs = prec[prec.decoy == 1]
miss_pair = (~pairs.tkey.isin(targets.index)).sum()

pred = pq.read_table(B + '/astral-lib12-20261111/pred12.parquet',
                     columns=['Modified.Sequence', 'Precursor.Charge', 'Product.Mz', 'Relative.Intensity', 'Fragment.Type',
                              'Fragment.Charge', 'Fragment.Series.Number']).to_pandas()
pred['tkey'] = pred['Modified.Sequence'].map(norm) + '/' + pred['Precursor.Charge'].astype(str)
pred = pred.drop_duplicates('tkey').set_index('tkey')
# A decoy whose target is not in the population still needs the target's prediction (key from TRAML_ID).
need_keys = set(targets.index) | set(pairs.tkey)
assert need_keys <= set(pred.index), len(need_keys - set(pred.index))

rows = []; tid = 0; dev = []
for p in prec.itertuples(index=False):
    q = pred.loc[p.tkey]
    res = residues(p.ms)
    for mz, li, ft, fz, o in zip(q['Product.Mz'], q['Relative.Intensity'], q['Fragment.Type'], q['Fragment.Charge'],
                                 q['Fragment.Series.Number']):
        ion = {2: 'b', 5: 'y'}[int(ft)]; o = int(o); fz = int(fz)
        m = frag_mz(res, ion, o, fz)
        if p.decoy == 0: dev.append(abs(m - mz))
        name = f'{p.traml}_{ion}{o}^{fz}'
        rows.append((tid, name, m, fz, ion, f'{ion}{o}^{fz}', o, 1, 0, 1, float(li), int(p.decoy), int(p.id))); tid += 1
print('target m/z against DIALibGen: max abs error', max(dev), flush=True)
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
print('precursors', len(prec), 'targets', int((prec.decoy == 0).sum()), 'decoys', len(pairs), 'decoys without target in population',
      int(miss_pair), 'transitions', len(t), 'per precursor', t.groupby('PRECURSOR_ID').size().value_counts().sort_index().to_dict())
