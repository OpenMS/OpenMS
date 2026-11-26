"""Label-free pair-calibration helpers. PROTEIN and biological annotations are never read here."""
import sqlite3
import numpy as np
import pandas as pd


def metadata(library):
    with sqlite3.connect(f'file:{library}?mode=ro', uri=True) as c:
        p = pd.read_sql('SELECT p.ID precursor_id, p.TRAML_ID name, p.DECOY decoy, p.CHARGE charge, p.PRECURSOR_MZ mz, p.LIBRARY_RT rt, p.LIBRARY_DRIFT_TIME im, pe.UNMODIFIED_SEQUENCE seq FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID JOIN PEPTIDE pe ON pe.ID=m.PEPTIDE_ID', c)
    assert not p.precursor_id.duplicated().any()
    p['key'] = np.where(p.decoy == 1, p.name.str.replace('^DECOY_', '', regex=True), p.name)
    return p


def pair_winners(table, meta):
    t = pd.read_csv(table, sep='\t', usecols=['precursor_id', 'score', 'is_decoy'])
    t = t[np.isfinite(t.score) & (t.score > -1e29)].sort_values(['score', 'precursor_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    t = t.merge(meta, on='precursor_id', validate='one_to_one')
    assert (t.is_decoy == t.decoy).all()
    t = t.sort_values(['score', 'decoy', 'precursor_id'], ascending=[False, False, True]).drop_duplicates('key')
    return t


def tdc_q(score, decoy):
    """Tie-grouped plus-one TDC, input and return in original row order."""
    score, decoy = np.asarray(score), np.asarray(decoy)
    assert np.isfinite(score).all()
    order = np.argsort(-score, kind='stable')
    s, d = score[order], decoy[order].astype(bool)
    if len(s) == 0:
        return np.empty(0)
    ends = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s)-1]
    dn, tn = np.cumsum(d)[ends], np.cumsum(~d)[ends]
    q = np.minimum(1, (dn + 1) / np.maximum(tn, 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    counts = np.diff(np.r_[-1, ends])
    result = np.empty(len(s))
    result[order] = np.repeat(q, counts)
    return result


def self_test():
    assert np.allclose(tdc_q([3,2,1], [0,0,0]), 1/3)
    assert np.allclose(tdc_q([3,2,1], [1,1,1]), 1)
    assert np.allclose(tdc_q([1,1,1], [0,1,0]), 1)
    a = tdc_q([4,3,2,1], [0,1,0,0])
    ix = np.array([3,1,0,2])
    assert np.allclose(a[ix], tdc_q(np.array([4,3,2,1])[ix], np.array([0,1,0,0])[ix]))
    print('PASS label-free TDC helper tests')


if __name__ == '__main__':
    self_test()
