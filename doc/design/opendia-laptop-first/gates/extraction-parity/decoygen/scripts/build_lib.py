"""PLAN_DECOYGEN step 3: one library per arm (A, B, C, D). Label-free; uses the library, the decoy map and the predictions.

Targets: the PAIRLIB target rows (Decoy = 0), copied verbatim.
Decoys: one per target precursor (TransitionGroupId T), named as the PAIRLIB decoys (TransitionGroupId DECOY_<T>, ProteinId
DECOY_<target ProteinId>) so pairing and every downstream script are unchanged. Its sequence is the arm's decoy of the
target's modified peptide (decoy_map.tsv); everything else is DIALibGen 0.13.1's prediction for that sequence at the
target's charge (decoys_p6.parquet): precursor m/z, the 6 fragments (m/z, annotation, relative intensity) and the
normalised RT. (The PAIRLIB targets' own NormalizedRetentionTime equals DIALibGen's RT; probe/ shows NCE 23 reproduces
the targets' 6 transitions best: median cosine 0.999.)
A target gets no decoy when its peptide has none in the arm (equals a target, or duplicates another decoy), when the
prediction lacks the decoy at that charge, or when the decoy precursor m/z is outside the library's 400-1000 window.
Also written: decoy12.tsv, the decoy's own 12-fragment prediction (decoys_p12.parquet, the pred12 settings), used by
build_pqp12_dg.py for the 12-fragment export library.
Usage: build_lib.py <arm>. Output: decoygen-20261113/arm<arm>/home/lib/library.tsv, decoy12.tsv, build_lib.txt (arm<arm>/home is
the arm's run root, bound onto the PAIRLIB path when the PAIRLIB chain runs)."""
import os, sys, numpy as np, pandas as pd, pyarrow.parquet as pq
ARM = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/decoygen-20261113'
O = f'{D}/arm{ARM}/home/lib'
os.makedirs(O, exist_ok=True)
log = open(O + '/build_lib.txt', 'w')
pr = lambda *a: (print(*a, flush=True), print(*a, file=log, flush=True))
norm = lambda s: s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')
ION = {2: 'b', 5: 'y'}

lib = pd.read_csv(B + '/pairlib-20261113/lib/pairlib.decoys.tsv', sep='\t', dtype=str, keep_default_na=False)
cols = list(lib.columns)
tgt = lib[lib.Decoy == '0']
pr('target transitions', len(tgt), 'target precursors', tgt.TransitionGroupId.nunique())
prec = tgt.drop_duplicates('TransitionGroupId')[['TransitionGroupId', 'ModifiedPeptideSequence', 'PrecursorCharge', 'ProteinId']]

dm = pd.read_csv(D + '/lib/decoy_map.tsv', sep='\t')
dm = dm[dm.arm == ARM].set_index('target_ms')
prec = prec.assign(decoy_ms=prec.ModifiedPeptideSequence.map(dm.decoy_ms), decoy_seq=prec.ModifiedPeptideSequence.map(dm.decoy_seq))
pr('precursors without a decoy peptide in this arm', int(prec.decoy_ms.isna().sum()))
prec = prec.dropna(subset=['decoy_ms'])


def load(path):
    p = pq.read_table(path).to_pandas()
    p['key'] = p['Modified.Sequence'].map(norm) + '/' + p['Precursor.Charge'].astype(str)
    return p.drop_duplicates('key').set_index('key')


p6 = load(D + '/lib/decoys_p6.parquet')
p12 = load(D + '/lib/decoys_p12.parquet')
prec['key'] = prec.decoy_ms + '/' + prec.PrecursorCharge
ok = prec.key.isin(p6.index)
pr('precursors whose decoy is missing from the 6-fragment prediction', int((~ok).sum()))
prec = prec[ok].copy()
six = p6.loc[prec.key, 'Product.Mz'].map(len).values == 6
pr('decoy precursors with fewer than 6 predicted fragments (dropped)', int((~six).sum()))
prec = prec[six].copy()
prec['pmz'] = p6.loc[prec.key, 'Precursor.Mz'].values
inwin = (prec.pmz >= 400) & (prec.pmz <= 1000)
pr('decoy precursors outside 400-1000 m/z (dropped)', int((~inwin).sum()))
prec = prec[inwin]
miss12 = (~prec.key.isin(p12.index)).sum()
pr('decoys missing from the 12-fragment prediction', int(miss12))

rows, rows12 = [], []
blank = {c: '' for c in cols}
for r in prec.itertuples(index=False):
    q = p6.loc[r.key]
    tg = 'DECOY_' + r.TransitionGroupId
    base = dict(blank, PrecursorMz=repr(float(q['Precursor.Mz'])), PrecursorCharge=r.PrecursorCharge,
                NormalizedRetentionTime=repr(float(q['RT'])), PeptideSequence=r.decoy_seq, ModifiedPeptideSequence=r.decoy_ms,
                PeptideGroupLabel='NA', ProteinId='DECOY_' + r.ProteinId, GeneName='NA', CollisionEnergy='-1.0',
                PrecursorIonMobility='-1.0', TransitionGroupId=tg, Decoy='1', DetectingTransition='1', IdentifyingTransition='0',
                QuantifyingTransition='1')
    it = np.asarray(q['Relative.Intensity'], float); it = it / it.max()
    for mz, li, ft, fz, o in zip(q['Product.Mz'], it, q['Fragment.Type'], q['Fragment.Charge'], q['Fragment.Series.Number']):
        ann = f'{ION[int(ft)]}{int(o)}^{int(fz)}'
        rows.append(dict(base, ProductMz=repr(float(mz)), ProductCharge=str(int(fz)), LibraryIntensity=repr(float(li)),
                         FragmentType=ION[int(ft)], FragmentSeriesNumber=str(int(o)), Annotation=ann, TransitionId=f'{tg}_{ann}'))
    if r.key in p12.index:
        q = p12.loc[r.key]
        for k, (mz, li, ft, fz, o) in enumerate(zip(q['Product.Mz'], q['Relative.Intensity'], q['Fragment.Type'], q['Fragment.Charge'],
                                                    q['Fragment.Series.Number'])):
            rows12.append((tg, r.decoy_ms, k, ION[int(ft)], int(o), int(fz), float(mz), float(li)))
dec = pd.DataFrame(rows, columns=cols)
pr('decoy precursors', len(prec), 'decoy transitions', len(dec), 'targets with a decoy %.4f' % (len(prec) / tgt.TransitionGroupId.nunique()))
pd.concat([tgt, dec]).to_csv(O + '/library.tsv', sep='\t', index=False)
pd.DataFrame(rows12, columns=['traml', 'decoy_ms', 'rank', 'ion', 'ord', 'fz', 'mz', 'li']).to_csv(O + '/decoy12.tsv', sep='\t', index=False)
# Library-only descriptors (no run data): decoy minus target precursor m/z, and the decoy RT against its target's.
t_pmz = tgt.drop_duplicates('TransitionGroupId').set_index('TransitionGroupId')
d = prec.set_index('TransitionGroupId')
dmz = d.pmz.values - t_pmz.loc[d.index, 'PrecursorMz'].astype(float).values
drt = p6.loc[d.key, 'RT'].values - t_pmz.loc[d.index, 'NormalizedRetentionTime'].astype(float).values
pr('decoy - target precursor m/z: |d| < 0.01 share %.3f, median |d| %.3f' % (np.mean(np.abs(dmz) < 0.01), np.median(np.abs(dmz))))
pr('decoy - target normalised RT: median |d| %.4f, RT range targets %.3f..%.3f decoys %.3f..%.3f' % (
    np.median(np.abs(drt)), t_pmz.NormalizedRetentionTime.astype(float).min(), t_pmz.NormalizedRetentionTime.astype(float).max(),
    p6.loc[d.key, 'RT'].min(), p6.loc[d.key, 'RT'].max()))
