"""Post-hoc control for block I (descriptive, not part of the PLAN_GAP ruling): is I's gain an artefact of the FDRBench
entrapment construction? Each entrapment peptide is a same-composition shuffle of a human peptide, so it is co-isolated
with its parent and can share fragments with it. If I mainly lowers entrapment scores through their parent claimant,
FDP_c would look better for a reason that does not exist in a real sample.
Variant J: features_I.py unchanged except that claimants with the same amino-acid composition and charge as the window's
precursor are ignored (label-free: composition comes from the library sequences, applied to targets and decoys alike).
Output featJ.npy (same 4 columns). Usage: features_J.py <raw> <root> <library12_pop.pqp> <iso_windows.tsv> <search.features.tsv>."""
import os, sqlite3, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(f'{HERE}/features_I.py').read()
hook = "    I, J = idx[I], idx[J]; keep = pid[I] != pid[J]; I, J = I[keep], J[keep]\n"
assert hook in src
src = src.replace(hook, "    I, J = idx[I], idx[J]; keep = (pid[I] != pid[J]) & (COMP[inv[I]] != COMP[inv[J]]); I, J = I[keep], J[keep]\n")
anchor = "mz = pmz.reindex(pop.precursor_id).values\n"
assert anchor in src
src = src.replace(anchor, anchor + """_seq = pd.read_sql('''select p.ID pid, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m
  on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID''', c).drop_duplicates('pid').set_index('pid')
_key = (_seq.seq.map(lambda s: ''.join(sorted(s))) + '/' + _seq.z.astype(str))
_codes = pd.Series(pd.factorize(_key)[0], index=_key.index)
COMP = _codes.reindex(uniq).fillna(-1).astype(int).values
COMP = np.where(COMP < 0, -1 - np.arange(len(COMP)), COMP)  # unknown composition: unique code, never equal
""")
src = src.replace("np.save(f'{O}/{RAW}/featI.npy', out)", "np.save(f'{O}/{RAW}/featJ.npy', out)").replace("RAW, 'featI'", "RAW, 'featJ'")
exec(compile(src, 'features_J', 'exec'))
